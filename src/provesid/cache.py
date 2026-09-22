"""
Advanced cache management for PROVESID package.

This module provides persistent, unlimited caching with size monitoring,
import/export capabilities, and automatic cache management.

Entries live under :func:`provesid.utils.user_cache_path` --- a real per-user
cache directory, not the system temp directory --- so a cache written today is
still there tomorrow. Set ``PROVESID_CACHE_DIR`` to put it somewhere else.

Three version numbers guard against serving an entry whose *shape* no longer
matches what the caller expects; see :data:`CACHE_KEY_VERSION`.
"""

import json
import pickle
import warnings
import hashlib
from pathlib import Path
from functools import wraps
from typing import Any, Dict, Optional, Union, Callable
from datetime import datetime

from .utils import user_cache_path


#: Key component bumped when something changes globally about what is stored or
#: how keys are derived, making *every* existing entry unreachable at once.
#: Two finer-grained versions sit beneath it, and the right one to bump is the
#: narrowest that covers the change:
#:
#: * a client's ``CACHE_SCHEMA_VERSION``, returned from its ``__cache_key__``,
#:   retires the entries of that one client (see
#:   :attr:`provesid.PubChemView.CACHE_SCHEMA_VERSION`);
#: * ``@cached(version=N)`` retires the entries of one function.
#:
#: A schema change is not something the cache can detect by itself: pickle
#: stores an instance's ``__dict__``, so an entry written before a dataclass
#: gained a field restores as an object missing that attribute, and every
#: caller reading it raises on a machine where only the package version
#: changed. Bumping a version is what makes those entries unreachable instead.
CACHE_KEY_VERSION = 1

#: The services with their own cache directory. This list is data: the
#: module-level cache functions all take ``service=`` rather than each service
#: owning a hand-written pair of functions.
CACHE_SERVICES = (
    'pubchem',
    'cas',
    'nci',
    'pubchemview',
    'classyfire',
    'opsin',
    'chebifier',
)


class _Miss:
    """Singleton marker for "no cache entry", distinct from a cached ``None``."""

    def __repr__(self) -> str:
        return "<cache miss>"


_MISS = _Miss()


def stable_key_part(obj: Any) -> Any:
    """
    Reduce a function argument to a JSON-serialisable, process-stable form.

    A cache key must come out identical for the same logical call in every
    process and for every client instance, otherwise a persistent entry written
    by one run is unreachable from the next. Python's default ``str()`` of an
    object embeds its memory address --- ``<PubChemView object at 0x7f...>`` ---
    so objects are reduced here to either the value they declare through
    ``__cache_key__()`` or their fully qualified class name.

    Args:
        obj: Any positional or keyword argument of a cached call. For a bound
            method this includes ``self``.

    Returns:
        A structure built only from strings, numbers, booleans, ``None``, lists
        and dicts, safe to pass to :func:`json.dumps` without ``default=str``.

    Note:
        Class-name reduction is right for a client instance, whose identity is
        its configuration, but it would make two *value* objects of the same
        class share a key. Any value-like object passed to a cached function
        must therefore declare ``__cache_key__``. No current caller passes one;
        cached functions take identifiers, strings and sequences.

    Examples:
        >>> stable_key_part({'cid': 2244, 'props': ('MolecularWeight',)})
        {'cid': 2244, 'props': ['MolecularWeight']}
    """
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, (list, tuple)):
        return [stable_key_part(item) for item in obj]
    if isinstance(obj, (set, frozenset)):
        return sorted((stable_key_part(item) for item in obj), key=repr)
    if isinstance(obj, dict):
        return {
            str(key): stable_key_part(value)
            for key, value in sorted(obj.items(), key=lambda kv: str(kv[0]))
        }

    declared = getattr(obj, '__cache_key__', None)
    if callable(declared):
        return stable_key_part(declared())
    return f"{type(obj).__module__}.{type(obj).__qualname__}"


def is_failure_result(result: Any) -> bool:
    """
    Report whether a return value is an in-band failure report.

    Several PROVESID clients signal an error by returning a dict with
    ``success: False`` and an ``error`` message rather than by raising. Storing
    one of those would turn a transient network failure into a permanent
    "no data" answer, so :func:`cached` never writes a value this function
    accepts.

    Args:
        result: The value returned by a cached function.

    Returns:
        True when ``result`` is a dict whose ``success`` key is ``False``.

    Examples:
        >>> is_failure_result({'success': False, 'error': 'Server busy'})
        True
        >>> is_failure_result({'success': True, 'CID': 2244})
        False
    """
    return isinstance(result, dict) and result.get('success') is False


def is_empty_result(result: Any) -> bool:
    """
    Report whether a lookup returned nothing.

    Absence is never cached. A service answers a genuinely empty lookup and a
    momentary refusal the same way often enough that the two cannot be told
    apart with confidence, and a wrongly cached "no data" is permanent whereas
    re-fetching a genuinely empty result costs one cheap request. Pass this to
    :func:`cached` as ``skip_if`` for any method that reports absence with an
    empty list, dict or DataFrame rather than by raising.

    Args:
        result: The value returned by a cached function.

    Returns:
        True when ``result`` is None or carries no items. A value with no
        length --- an int, a dataclass --- is a real result.

    Examples:
        >>> is_empty_result([])
        True
        >>> is_empty_result(None)
        True
        >>> is_empty_result(['3.47'])
        False
    """
    if result is None:
        return True
    try:
        return len(result) == 0
    except TypeError:
        # Not a sized value; treat it as a real result.
        return False


class CacheManager:
    """
    Advanced cache manager with persistent storage, size monitoring, and import/export.

    Features:
    - Unlimited cache by default
    - Persistent storage across sessions, in a real cache directory rather than
      the system temp directory
    - Size monitoring with warnings at 5GB
    - Import/export functionality
    - Cache statistics and management

    Examples:
        >>> import tempfile
        >>> with tempfile.TemporaryDirectory() as tmp:
        ...     manager = CacheManager(cache_dir=tmp)
        ...     manager.set('demo', (2244,), {}, {'cid': 2244})
        ...     manager.get('demo', (2244,), {})
        (True, {'cid': 2244})
    """

    def __init__(self,
                 cache_dir: Optional[str] = None,
                 service_name: Optional[str] = None,
                 max_size_gb: float = 5.0,
                 enable_warnings: bool = True):
        """
        Initialize the cache manager.

        Args:
            cache_dir: Directory to store cache files. When None, a per-user
                cache directory from :func:`provesid.utils.user_cache_path` is
                used, with ``service_name`` as a subdirectory if given. The
                system temp directory is deliberately not the default: ``/tmp``
                is cleared on boot on most Linux systems, which made the
                "persistent" cache last only until the next reboot.
            service_name: Name of the service for service-specific caching
                (one of :data:`CACHE_SERVICES`).
            max_size_gb: Size threshold in GB for warnings (default: 5.0)
            enable_warnings: Whether to enable size warnings (default: True)
        """
        self.service_name = service_name
        self.max_size_gb = max_size_gb
        self.enable_warnings = enable_warnings
        self._last_size_check = None
        self._size_check_interval = 100  # Check size every 100 cache operations
        self._operation_count = 0

        # Set up cache directory
        if cache_dir is None:
            cache_dir = (user_cache_path(service_name, ensure_exists=False)
                         if service_name else
                         user_cache_path(ensure_exists=False))

        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Cache storage - in-memory cache backed by persistent storage
        self._memory_cache: Dict[str, Any] = {}
        self._cache_metadata: Dict[str, Dict] = {}

        # Load existing cache metadata
        self._load_metadata()

    def _load_metadata(self):
        """Load cache metadata from persistent storage."""
        metadata_file = self.cache_dir / 'metadata.json'
        if metadata_file.exists():
            try:
                with open(metadata_file, 'r') as f:
                    self._cache_metadata = json.load(f)
            except Exception as e:
                warnings.warn(f"Could not load cache metadata: {e}")
                self._cache_metadata = {}

    def _save_metadata(self):
        """Save cache metadata to persistent storage."""
        metadata_file = self.cache_dir / 'metadata.json'
        try:
            with open(metadata_file, 'w') as f:
                json.dump(self._cache_metadata, f, indent=2, default=str)
        except Exception as e:
            warnings.warn(f"Could not save cache metadata: {e}")

    def _get_cache_key(self, func_name: str, args: tuple, kwargs: dict,
                       version: int = 1) -> str:
        """
        Generate a unique, process-stable cache key for a function call.

        Args:
            func_name: Fully qualified name of the cached function.
            args: Positional arguments of the call. For a bound method the first
                entry is the client instance; :func:`stable_key_part` reduces it
                to its ``__cache_key__`` or class name so the key does not
                depend on the instance's memory address.
            kwargs: Keyword arguments of the call.
            version: Shape version of what this one function returns, from
                ``@cached(version=...)``. Bumping it retires that function's
                entries and nothing else. :data:`CACHE_KEY_VERSION` is folded
                in as well, so a global bump retires everything.

        Returns:
            Hex SHA-256 digest of the normalised call signature.
        """
        key_data = {
            'cache_version': CACHE_KEY_VERSION,
            'function': func_name,
            'function_version': version,
            'args': stable_key_part(list(args)),
            'kwargs': stable_key_part(dict(kwargs)),
        }
        key_str = json.dumps(key_data, sort_keys=True)
        return hashlib.sha256(key_str.encode()).hexdigest()

    def _get_cache_file_path(self, cache_key: str) -> Path:
        """Get the file path for a cache entry."""
        return self.cache_dir / f"{cache_key}.pkl"

    def _load_from_disk(self, cache_key: str) -> Any:
        """
        Load a cache entry from disk.

        Args:
            cache_key: Key produced by :meth:`_get_cache_key`.

        Returns:
            The stored value, or the ``_MISS`` sentinel when there is no
            readable entry. A stored ``None`` comes back as ``None``, so a
            function that legitimately returns ``None`` still caches.
        """
        cache_file = self._get_cache_file_path(cache_key)
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                warnings.warn(f"Could not load cache entry {cache_key}: {e}")
        return _MISS

    def _save_to_disk(self, cache_key: str, value: Any):
        """Save a cache entry to disk."""
        cache_file = self._get_cache_file_path(cache_key)
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(value, f)

            # Update metadata
            self._cache_metadata[cache_key] = {
                'created': datetime.now().isoformat(),
                'size': cache_file.stat().st_size,
                'file': str(cache_file)
            }

            # Periodically save metadata and check size
            self._operation_count += 1
            if self._operation_count % self._size_check_interval == 0:
                self._save_metadata()
                self._check_cache_size()

        except Exception as e:
            warnings.warn(f"Could not save cache entry {cache_key}: {e}")

    def get(self, func_name: str, args: tuple, kwargs: dict,
            version: int = 1) -> tuple:
        """
        Get a cached result.

        Args:
            func_name: Fully qualified name of the cached function.
            args: Positional arguments of the call.
            kwargs: Keyword arguments of the call.
            version: Shape version of the function's return value; must match
                the one used when the entry was written, or it is a miss.

        Returns:
            (found: bool, value: Any) - Tuple indicating if value was found and the value.
            ``found`` is what tells a stored ``None`` from a miss.

        Examples:
            >>> import tempfile
            >>> manager = CacheManager(cache_dir=tempfile.mkdtemp())
            >>> manager.get('demo', (2244,), {})
            (False, None)
            >>> manager.set('demo', (2244,), {}, None)
            >>> manager.get('demo', (2244,), {})
            (True, None)
            >>> manager.get('demo', (2244,), {}, version=2)     # another shape
            (False, None)
        """
        cache_key = self._get_cache_key(func_name, args, kwargs, version)

        # Check memory cache first
        if cache_key in self._memory_cache:
            return True, self._memory_cache[cache_key]

        # Check disk cache
        value = self._load_from_disk(cache_key)
        if value is not _MISS:
            # Load into memory cache
            self._memory_cache[cache_key] = value
            return True, value

        return False, None

    def set(self, func_name: str, args: tuple, kwargs: dict, value: Any,
            version: int = 1):
        """
        Store a result under the key for this call.

        Args:
            func_name: Fully qualified name of the cached function.
            args: Positional arguments of the call.
            kwargs: Keyword arguments of the call.
            value: The value to store. It must pickle.
            version: Shape version of the function's return value.

        Note:
            A value that cannot be written to disk stays in memory for this
            process, with a warning.

        Examples:
            >>> import tempfile
            >>> manager = CacheManager(cache_dir=tempfile.mkdtemp())
            >>> manager.set('demo', ('aspirin',), {'exact': True}, [2244])
            >>> manager.get('demo', ('aspirin',), {'exact': True})
            (True, [2244])
            >>> manager.get('demo', ('aspirin',), {'exact': False})
            (False, None)
        """
        cache_key = self._get_cache_key(func_name, args, kwargs, version)

        # Store in memory cache
        self._memory_cache[cache_key] = value

        # Store on disk
        self._save_to_disk(cache_key, value)

    def clear(self):
        """
        Delete every entry of this cache, in memory and on disk.

        Only this manager's directory is touched; other services' caches are
        kept.

        Examples:
            >>> import tempfile
            >>> manager = CacheManager(cache_dir=tempfile.mkdtemp())
            >>> manager.set('demo', (1,), {}, 'one')
            >>> manager.clear()
            >>> manager.get('demo', (1,), {})
            (False, None)
        """
        # Clear memory cache
        self._memory_cache.clear()

        # Clear disk cache
        for cache_file in self.cache_dir.glob("*.pkl"):
            try:
                cache_file.unlink()
            except Exception as e:
                warnings.warn(f"Could not delete cache file {cache_file}: {e}")

        # Clear metadata
        self._cache_metadata.clear()
        self._save_metadata()

    def get_cache_size(self) -> Dict[str, Union[int, float]]:
        """
        Get cache size information.

        Returns:
            Dictionary with size information in bytes, MB, and GB: ``bytes``,
            ``mb``, ``gb`` and ``files``, counted from the entry files on disk.

        Examples:
            >>> import tempfile
            >>> manager = CacheManager(cache_dir=tempfile.mkdtemp())
            >>> manager.set('demo', (1,), {}, 'one')
            >>> size = manager.get_cache_size()
            >>> size['files'], size['bytes'] > 0
            (1, True)
        """
        total_size = 0
        file_count = 0

        for cache_file in self.cache_dir.glob("*.pkl"):
            try:
                total_size += cache_file.stat().st_size
                file_count += 1
            except Exception:
                continue

        return {
            'bytes': total_size,
            'mb': total_size / (1024 * 1024),
            'gb': total_size / (1024 * 1024 * 1024),
            'files': file_count
        }

    def _check_cache_size(self):
        """Check cache size and issue warnings if necessary."""
        if not self.enable_warnings:
            return

        size_info = self.get_cache_size()
        if size_info['gb'] > self.max_size_gb:
            warnings.warn(
                f"Cache size ({size_info['gb']:.2f} GB) exceeds warning threshold "
                f"({self.max_size_gb} GB). Consider using export_cache() to backup "
                f"and clear() to free space, or increase the threshold.",
                UserWarning
            )

    def get_cache_info(self) -> Dict[str, Any]:
        """
        Get comprehensive cache information.

        Returns:
            dict: ``cache_directory``; ``memory_entries``, the entries this
            process holds in memory; ``disk_entries`` and ``file_count``, the
            entries on disk, whichever process wrote them; ``total_size_bytes``,
            ``total_size_mb`` and ``total_size_gb``; ``warning_threshold_gb``
            and ``warnings_enabled``.

        Examples:
            >>> import tempfile
            >>> manager = CacheManager(cache_dir=tempfile.mkdtemp())
            >>> manager.set('demo', (1,), {}, 'one')
            >>> info = manager.get_cache_info()
            >>> info['memory_entries'], info['disk_entries'], info['warning_threshold_gb']
            (1, 1, 5.0)
        """
        size_info = self.get_cache_size()

        return {
            'cache_directory': str(self.cache_dir),
            'memory_entries': len(self._memory_cache),
            # The files, not the metadata: metadata.json is flushed every
            # hundred writes, so a short-lived process never records its own.
            'disk_entries': size_info['files'],
            'total_size_bytes': size_info['bytes'],
            'total_size_mb': size_info['mb'],
            'total_size_gb': size_info['gb'],
            'file_count': size_info['files'],
            'warning_threshold_gb': self.max_size_gb,
            'warnings_enabled': self.enable_warnings
        }

    def export_cache(self, export_path: str, format: str = 'pickle') -> bool:
        """
        Export cache data to a file.

        Every entry on disk is exported, whichever process wrote it.

        Args:
            export_path: Path to export file
            format: Export format ('pickle', 'json'). JSON writes values that
                are not JSON types as their ``str()``, so only a pickle export
                imports back unchanged.

        Returns:
            True if export successful, False otherwise (with a warning)

        Examples:
            >>> import os, tempfile
            >>> manager = CacheManager(cache_dir=tempfile.mkdtemp())
            >>> manager.set('demo', (1,), {}, 'one')
            >>> backup = os.path.join(tempfile.mkdtemp(), 'backup.pkl')
            >>> manager.export_cache(backup)
            True
        """
        try:
            export_data = {
                'metadata': self._cache_metadata,
                'cache_data': {},
                'export_info': {
                    'timestamp': datetime.now().isoformat(),
                    'version': '1.0',
                    'format': format
                }
            }

            # Every entry file, not just those in the metadata, which is
            # flushed only every hundred writes and so misses whatever a
            # short-lived process stored.
            for cache_file in sorted(self.cache_dir.glob("*.pkl")):
                cache_key = cache_file.stem
                value = self._load_from_disk(cache_key)
                # `is not _MISS`, not `is not None`: a function may legitimately
                # return None, and that entry has to survive the round trip.
                if value is not _MISS:
                    export_data['cache_data'][cache_key] = value

            if format == 'pickle':
                with open(export_path, 'wb') as f:
                    pickle.dump(export_data, f)
            elif format == 'json':
                with open(export_path, 'w') as f:
                    json.dump(export_data, f, indent=2, default=str)
            else:
                raise ValueError(f"Unsupported format: {format}")

            return True

        except Exception as e:
            warnings.warn(f"Cache export failed: {e}")
            return False

    def import_cache(self, import_path: str, merge: bool = True) -> bool:
        """
        Import cache data from a file.

        Imported entries are written to disk and found by the next
        :meth:`get`; the keys are kept as exported, so an export imports
        into any cache whose functions and versions match.

        Args:
            import_path: Path to import file. A name ending ``.json`` is read
                as JSON, anything else as pickle.
            merge: If True, merge with existing cache. If False, replace existing cache.

        Returns:
            True if import successful, False otherwise (with a warning)

        Examples:
            >>> import os, tempfile
            >>> source = CacheManager(cache_dir=tempfile.mkdtemp())
            >>> source.set('demo', (1,), {}, 'one')
            >>> backup = os.path.join(tempfile.mkdtemp(), 'backup.pkl')
            >>> source.export_cache(backup)
            True
            >>> target = CacheManager(cache_dir=tempfile.mkdtemp())
            >>> target.import_cache(backup)
            True
            >>> target.get('demo', (1,), {})
            (True, 'one')
        """
        try:
            # Determine format from file extension or content
            if import_path.endswith('.json'):
                with open(import_path, 'r') as f:
                    import_data = json.load(f)
            else:
                with open(import_path, 'rb') as f:
                    import_data = pickle.load(f)

            if not merge:
                self.clear()

            # Import cache data
            metadata = import_data.get('metadata', {})
            cache_data = import_data.get('cache_data', {})

            for cache_key, value in cache_data.items():
                # Save to disk
                cache_file = self._get_cache_file_path(cache_key)
                with open(cache_file, 'wb') as f:
                    pickle.dump(value, f)

                # Update metadata
                if cache_key in metadata:
                    self._cache_metadata[cache_key] = metadata[cache_key]
                else:
                    self._cache_metadata[cache_key] = {
                        'created': datetime.now().isoformat(),
                        'size': cache_file.stat().st_size,
                        'file': str(cache_file),
                        'imported': True
                    }

            self._save_metadata()
            return True

        except Exception as e:
            warnings.warn(f"Cache import failed: {e}")
            return False

# Cache managers are built on first use, not on import: constructing one
# creates its directory, and importing provesid should not leave eight empty
# directories behind on a machine that never calls a single service.
_global_cache: Optional[CacheManager] = None
_service_caches: Dict[str, CacheManager] = {}


def get_service_cache(service: Optional[str] = None) -> CacheManager:
    """
    Return the cache manager for one service, building it on first use.

    Args:
        service: One of :data:`CACHE_SERVICES`, or None for the global cache
            that holds everything not attributed to a service.

    Returns:
        The shared :class:`CacheManager` for that service. The same object is
        returned on every call, so two clients of the same service share
        entries and share the in-memory tier.

    Raises:
        ValueError: If ``service`` is neither in :data:`CACHE_SERVICES` nor
            already registered. Silently falling back to the global cache
            would scatter a typo's entries somewhere the matching
            ``clear_cache`` call never looks. A name registered directly in
            ``_service_caches`` --- which is how tests inject a throwaway cache
            --- counts as known.

    Examples:
        >>> get_service_cache('pubchem') is get_service_cache('pubchem')
        True
    """
    global _global_cache

    if service is None:
        if _global_cache is None:
            _global_cache = CacheManager()
        return _global_cache

    if service not in _service_caches:
        _require_known_service(service)
        _service_caches[service] = CacheManager(service_name=service)
    return _service_caches[service]


def _require_known_service(service: str) -> None:
    """Raise unless ``service`` is one this module is prepared to cache for."""
    if service not in CACHE_SERVICES and service not in _service_caches:
        raise ValueError(
            f"Unknown cache service: {service!r}. "
            f"Known services: {', '.join(CACHE_SERVICES)}"
        )


def cached(func: Callable = None, *, service: Optional[str] = None,
           skip_if: Optional[Callable[[Any], bool]] = None,
           version: int = 1) -> Callable:
    """
    Decorator for unlimited caching with persistent storage.

    This replaces the standard ``@lru_cache`` decorator with unlimited caching,
    persistent storage, and size monitoring. Keys are process-stable, so an
    entry written by one run is found by the next (see
    :func:`stable_key_part`).

    Failed lookups are never written. An exception propagates uncached, and a
    return value that :func:`is_failure_result` recognises as an in-band failure
    report is returned to the caller but not stored --- a transient HTTP 429,
    503 or timeout must not become a permanent "no data" answer.

    Args:
        func: The function to cache, when used as a bare ``@cached``.
        service: Optional service name for service-specific caching; one of
            :data:`CACHE_SERVICES`. Omit it to use the global cache.
        skip_if: Optional predicate applied to the return value; when it returns
            True the value is handed back but not stored. Use it for functions
            that report a failed lookup with something other than a
            ``success: False`` dict.
        version: Shape version of this function's return value, part of every
            key it writes. Bump it when the function starts returning a
            different *shape* --- a new dataclass field, a dict key renamed,
            a list of strings becoming a list of objects --- so that entries of
            the old shape become unreachable instead of being unpickled into
            something the caller cannot read. Bump the client's
            ``CACHE_SCHEMA_VERSION`` instead when the change affects every
            method of one client, and :data:`CACHE_KEY_VERSION` when it affects
            everything.

    Returns:
        The decorated function, carrying ``cache_clear`` and ``cache_info``
        attributes bound to the selected cache.

    Raises:
        ValueError: At decoration time, if ``service`` is not a known service.

    Note:
        ``use_cache=False`` --- either as a constructor flag on the client or as
        a keyword on the call --- means "do not *read* the cache"; a successful
        result is still written so later calls benefit.

    Examples:
        >>> calls = []
        >>> @cached(service='pubchem', version=2)
        ... def fetch(cid):
        ...     calls.append(cid)
        ...     return {'success': True, 'cid': cid}
        >>> fetch.cache_clear()
        >>> fetch(2244)['cid'], fetch(2244)['cid'], calls
        (2244, 2244, [2244])
    """
    def decorator(f: Callable) -> Callable:
        @wraps(f)
        def wrapper(*args, **kwargs):
            func_name = f"{f.__module__}.{f.__qualname__}"
            cache_manager = get_service_cache(service)

            # Check if we should use cache
            use_cache = True

            # For instance methods, check self.use_cache
            if args and hasattr(args[0], 'use_cache'):
                use_cache = getattr(args[0], 'use_cache', True)

            # For static/class methods with use_cache parameter in kwargs
            if 'use_cache' in kwargs:
                use_cache = kwargs.pop('use_cache')  # Remove from kwargs to avoid passing to function

            # Try to get from cache only if use_cache is True
            if use_cache:
                found, value = cache_manager.get(func_name, args, kwargs, version)
                if found:
                    return value

            # Call the function. A raised exception propagates and is not cached.
            result = f(*args, **kwargs)

            # Never persist a failed lookup, even when use_cache is False:
            # caching a transient network error would make it permanent.
            if is_failure_result(result) or (skip_if is not None and skip_if(result)):
                return result

            cache_manager.set(func_name, args, kwargs, result, version)

            return result

        # Reject an unknown service now rather than on the first call, and
        # bind the management helpers lazily so decoration creates nothing.
        if service is not None:
            _require_known_service(service)
        wrapper.cache_clear = lambda: get_service_cache(service).clear()
        wrapper.cache_info = lambda: get_service_cache(service).get_cache_info()

        return wrapper

    # Support both @cached and @cached(service='name') syntax
    if func is None:
        return decorator
    else:
        return decorator(func)


# ---------------------------------------------------------------------------
# Cache management. Every function takes ``service=`` rather than each service
# owning a hand-written twin: the service list is data (:data:`CACHE_SERVICES`),
# and fourteen near-identical one-line functions were fourteen places to forget
# when an eighth service arrives.
# ---------------------------------------------------------------------------

def clear_cache(service: Optional[str] = None, all_services: bool = False):
    """
    Delete every cached entry for one service, or for the global cache.

    Args:
        service: One of :data:`CACHE_SERVICES`, or None for the global cache.
        all_services: When True, clear every service cache as well as the
            global one, and ignore ``service``.

    Raises:
        ValueError: If ``service`` names no known service.

    Examples:
        >>> clear_cache(service='pubchem')
    """
    if all_services:
        get_service_cache(None).clear()
        for name in CACHE_SERVICES:
            get_service_cache(name).clear()
        return
    get_service_cache(service).clear()


def get_cache_info(service: Optional[str] = None) -> Dict[str, Any]:
    """
    Report directory, entry counts and size for one cache.

    Args:
        service: One of :data:`CACHE_SERVICES`, or None for the global cache.

    Returns:
        Dict with ``cache_directory``, ``memory_entries``, ``disk_entries``,
        ``total_size_bytes``/``_mb``/``_gb``, ``file_count``,
        ``warning_threshold_gb`` and ``warnings_enabled``.

    Raises:
        ValueError: If ``service`` names no known service.

    Examples:
        >>> sorted(get_cache_info(service='opsin'))[:2]
        ['cache_directory', 'disk_entries']
    """
    return get_service_cache(service).get_cache_info()


def get_all_cache_info() -> Dict[str, Dict[str, Any]]:
    """
    Report :func:`get_cache_info` for the global cache and every service.

    Returns:
        Dict keyed by ``'global'`` and by each name in :data:`CACHE_SERVICES`.

    Examples:
        >>> 'pubchem' in get_all_cache_info()
        True
    """
    info = {'global': get_service_cache(None).get_cache_info()}
    for name in CACHE_SERVICES:
        info[name] = get_service_cache(name).get_cache_info()
    return info


def export_cache(export_path: str, format: str = 'pickle',
                 service: Optional[str] = None) -> bool:
    """
    Write one cache's entries to a file so they can be shared or archived.

    Args:
        export_path: Destination file.
        format: ``'pickle'`` (default, preserves Python objects) or ``'json'``.
        service: One of :data:`CACHE_SERVICES`, or None for the global cache.

    Returns:
        True on success, False if the export failed (a warning says why).

    Raises:
        ValueError: If ``service`` names no known service.

    Examples:
        >>> export_cache('cache.pkl', service='pubchem')    # doctest: +SKIP
        True
    """
    return get_service_cache(service).export_cache(export_path, format)


def import_cache(import_path: str, merge: bool = True,
                 service: Optional[str] = None) -> bool:
    """
    Load entries exported by :func:`export_cache` into one cache.

    Args:
        import_path: File written by :func:`export_cache`.
        merge: When True (default), keep existing entries; when False, clear
            the cache first.
        service: One of :data:`CACHE_SERVICES`, or None for the global cache.

    Returns:
        True on success, False if the import failed (a warning says why).

    Raises:
        ValueError: If ``service`` names no known service.

    Note:
        Entries carry the key they were written under, so an export made before
        a version bump imports cleanly and is simply never read.

    Examples:
        >>> import_cache('cache.pkl', service='pubchem')    # doctest: +SKIP
        True
    """
    return get_service_cache(service).import_cache(import_path, merge)


def get_cache_size(service: Optional[str] = None) -> Dict[str, Union[int, float]]:
    """
    Measure one cache on disk.

    Args:
        service: One of :data:`CACHE_SERVICES`, or None for the global cache.

    Returns:
        Dict with ``bytes``, ``mb``, ``gb`` and ``files``.

    Raises:
        ValueError: If ``service`` names no known service.

    Examples:
        >>> sorted(get_cache_size())
        ['bytes', 'files', 'gb', 'mb']
    """
    return get_service_cache(service).get_cache_size()


def set_cache_warning_threshold(size_gb: float, service: Optional[str] = None):
    """
    Set the size at which a cache starts warning, in GB.

    Args:
        size_gb: New threshold.
        service: One of :data:`CACHE_SERVICES`, or None for the global cache.
            The setting is per cache and is not remembered across processes.

    Raises:
        ValueError: If ``service`` names no known service.

    Examples:
        >>> set_cache_warning_threshold(10.0, service='pubchem')  # doctest: +SKIP
    """
    get_service_cache(service).max_size_gb = size_gb


def enable_cache_warnings(enabled: bool = True, service: Optional[str] = None):
    """
    Turn size warnings on or off for a cache.

    Args:
        enabled: True to warn past the threshold, False to stay quiet.
        service: One of :data:`CACHE_SERVICES`, or None for the global cache.

    Raises:
        ValueError: If ``service`` names no known service.

    Examples:
        >>> enable_cache_warnings(False, service='pubchem')  # doctest: +SKIP
    """
    get_service_cache(service).enable_warnings = enabled
