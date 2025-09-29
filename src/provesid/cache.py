"""
Advanced cache management for PROVESID package.

This module provides persistent, unlimited caching with size monitoring,
import/export capabilities, and automatic cache management.
"""

import os
import json
import pickle
import warnings
import hashlib
import sys
from pathlib import Path
from functools import wraps
from typing import Any, Dict, Optional, Union, Callable
from datetime import datetime
import tempfile

class CacheManager:
    """
    Advanced cache manager with persistent storage, size monitoring, and import/export.
    
    Features:
    - Unlimited cache by default
    - Persistent storage across sessions
    - Size monitoring with warnings at 5GB
    - Import/export functionality
    - Cache statistics and management
    """
    
    def __init__(self, 
                 cache_dir: Optional[str] = None,
                 max_size_gb: float = 5.0,
                 enable_warnings: bool = True):
        """
        Initialize the cache manager.
        
        Args:
            cache_dir: Directory to store cache files. If None, uses system temp directory.
            max_size_gb: Size threshold in GB for warnings (default: 5.0)
            enable_warnings: Whether to enable size warnings (default: True)
        """
        self.max_size_gb = max_size_gb
        self.enable_warnings = enable_warnings
        self._last_size_check = None
        self._size_check_interval = 100  # Check size every 100 cache operations
        self._operation_count = 0
        
        # Set up cache directory
        if cache_dir is None:
            cache_dir = os.path.join(tempfile.gettempdir(), 'provesid_cache')
        
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
    
    def _get_cache_key(self, func_name: str, args: tuple, kwargs: dict) -> str:
        """Generate a unique cache key for function call."""
        # Create deterministic hash from function name and arguments
        key_data = {
            'function': func_name,
            'args': args,
            'kwargs': sorted(kwargs.items()) if kwargs else {}
        }
        key_str = json.dumps(key_data, sort_keys=True, default=str)
        return hashlib.sha256(key_str.encode()).hexdigest()
    
    def _get_cache_file_path(self, cache_key: str) -> Path:
        """Get the file path for a cache entry."""
        return self.cache_dir / f"{cache_key}.pkl"
    
    def _load_from_disk(self, cache_key: str) -> Any:
        """Load a cache entry from disk."""
        cache_file = self._get_cache_file_path(cache_key)
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                warnings.warn(f"Could not load cache entry {cache_key}: {e}")
        return None
    
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
    
    def get(self, func_name: str, args: tuple, kwargs: dict) -> tuple:
        """
        Get a cached result.
        
        Returns:
            (found: bool, value: Any) - Tuple indicating if value was found and the value
        """
        cache_key = self._get_cache_key(func_name, args, kwargs)
        
        # Check memory cache first
        if cache_key in self._memory_cache:
            return True, self._memory_cache[cache_key]
        
        # Check disk cache
        value = self._load_from_disk(cache_key)
        if value is not None:
            # Load into memory cache
            self._memory_cache[cache_key] = value
            return True, value
        
        return False, None
    
    def set(self, func_name: str, args: tuple, kwargs: dict, value: Any):
        """Set a cached result."""
        cache_key = self._get_cache_key(func_name, args, kwargs)
        
        # Store in memory cache
        self._memory_cache[cache_key] = value
        
        # Store on disk
        self._save_to_disk(cache_key, value)
    
    def clear(self):
        """Clear all cached data."""
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
            Dictionary with size information in bytes, MB, and GB
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
        """Get comprehensive cache information."""
        size_info = self.get_cache_size()
        
        return {
            'cache_directory': str(self.cache_dir),
            'memory_entries': len(self._memory_cache),
            'disk_entries': len(self._cache_metadata),
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
        
        Args:
            export_path: Path to export file
            format: Export format ('pickle', 'json')
            
        Returns:
            True if export successful, False otherwise
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
            
            # Load all cache data
            for cache_key in self._cache_metadata:
                value = self._load_from_disk(cache_key)
                if value is not None:
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
        
        Args:
            import_path: Path to import file
            merge: If True, merge with existing cache. If False, replace existing cache.
            
        Returns:
            True if import successful, False otherwise
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


# Global cache manager instance
_global_cache = CacheManager()

def cached(func: Callable) -> Callable:
    """
    Decorator for unlimited caching with persistent storage.
    
    This replaces the standard @lru_cache decorator with unlimited caching,
    persistent storage, and size monitoring.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        func_name = f"{func.__module__}.{func.__qualname__}"
        
        # Try to get from cache
        found, value = _global_cache.get(func_name, args, kwargs)
        if found:
            return value
        
        # Call function and cache result
        result = func(*args, **kwargs)
        _global_cache.set(func_name, args, kwargs, result)
        
        return result
    
    # Add cache management methods to the decorated function
    wrapper.cache_clear = lambda: _global_cache.clear()
    wrapper.cache_info = lambda: _global_cache.get_cache_info()
    
    return wrapper

# Convenience functions for global cache management
def clear_cache():
    """Clear all cached data from the global cache."""
    _global_cache.clear()

def get_cache_info() -> Dict[str, Any]:
    """Get information about the global cache."""
    return _global_cache.get_cache_info()

def export_cache(export_path: str, format: str = 'pickle') -> bool:
    """Export global cache to a file."""
    return _global_cache.export_cache(export_path, format)

def import_cache(import_path: str, merge: bool = True) -> bool:
    """Import cache data from a file."""
    return _global_cache.import_cache(import_path, merge)

def get_cache_size() -> Dict[str, Union[int, float]]:
    """Get cache size information."""
    return _global_cache.get_cache_size()

def set_cache_warning_threshold(size_gb: float):
    """Set the cache size warning threshold in GB."""
    _global_cache.max_size_gb = size_gb

def enable_cache_warnings(enabled: bool = True):
    """Enable or disable cache size warnings."""
    _global_cache.enable_warnings = enabled