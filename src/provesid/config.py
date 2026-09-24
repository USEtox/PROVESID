"""
Configuration management for PROVESID API keys and settings.
Provides persistent storage for API keys and user preferences.

The keys live in ``config.json`` under ``$XDG_CONFIG_HOME/provesid`` (by
default ``~/.config/provesid``), or ``%APPDATA%\\PROVESID`` on Windows, as
plain text. Only CAS Common Chemistry needs one today:
[`CASCommonChem`][provesid.cascommonchem.CASCommonChem] reads it when no key or key
file is passed, and before the ``CCC_API_KEY`` and ``CAS_API_KEY``
environment variables --- so a stored key wins over the environment.

Examples:
    >>> from provesid.config import get_config_manager
    >>> get_config_manager().config_file.name
    'config.json'
"""

import os
import json
import logging
from pathlib import Path
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)


class ConfigManager:
    """
    Manages persistent configuration for PROVESID.

    Every read goes to the file, so a key set in another process is seen at
    once. Use [`get_config_manager`][provesid.config.get_config_manager] for
    the shared instance.

    Examples:
        >>> manager = ConfigManager()
        >>> manager.set_api_key("example", "not-a-real-key")
        >>> manager.get_api_key("example")
        'not-a-real-key'
        >>> manager.remove_api_key("example")
        True
    """

    def __init__(self):
        """
        Initialize configuration manager with default paths.

        Creates the configuration directory if it is missing; a failure to
        create it is logged, not raised.

        Examples:
            >>> ConfigManager().config_dir.name in ("provesid", "PROVESID")
            True
        """
        self.config_dir = self._get_config_directory()
        self.config_file = self.config_dir / "config.json"
        self._ensure_config_directory()

    def _get_config_directory(self) -> Path:
        """Get the appropriate configuration directory for the current OS"""
        if os.name == 'nt':  # Windows
            config_root = Path(os.environ.get('APPDATA', os.path.expanduser('~'))) / 'PROVESID'
        else:  # Unix-like (Linux, macOS)
            config_root = Path(os.environ.get('XDG_CONFIG_HOME', os.path.expanduser('~/.config'))) / 'provesid'

        return config_root

    def _ensure_config_directory(self):
        """Create configuration directory if it doesn't exist"""
        try:
            self.config_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logging.warning(f"Could not create config directory {self.config_dir}: {e}")

    def load_config(self) -> Dict[str, Any]:
        """
        Load configuration from file.

        Returns:
            (dict): The file's contents; empty when there is no file or it cannot
            be read (the latter logged).

        Examples:
            >>> isinstance(ConfigManager().load_config(), dict)
            True
        """
        if not self.config_file.exists():
            return {}

        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            logging.warning(f"Error loading config from {self.config_file}: {e}")
            return {}

    def save_config(self, config: Dict[str, Any]):
        """
        Save configuration to file, replacing what is there.

        Args:
            config: The whole configuration. Keys not in it are lost; to change
                one entry, load, modify and save.

        Note:
            A failure to write is logged at ERROR, not raised.

        Examples:
            >>> manager = ConfigManager()
            >>> config = manager.load_config()
            >>> manager.save_config(config)
        """
        try:
            with open(self.config_file, 'w', encoding='utf-8') as f:
                json.dump(config, f, indent=2, ensure_ascii=False)
        except Exception as e:
            logging.error(f"Error saving config to {self.config_file}: {e}")

    def get_api_key(self, service: str) -> Optional[str]:
        """
        Get API key for a specific service.

        Args:
            service: Service name, e.g. ``"cas"``.

        Returns:
            The stored key, or None when there is none.

        Examples:
            >>> ConfigManager().get_api_key("no-such-service") is None
            True
        """
        config = self.load_config()
        api_keys = config.get('api_keys', {})
        return api_keys.get(service)

    def set_api_key(self, service: str, api_key: str):
        """
        Set API key for a specific service, replacing any stored one.

        Args:
            service: Service name, e.g. ``"cas"``.
            api_key: The key; surrounding whitespace is stripped.

        Examples:
            >>> manager = ConfigManager()
            >>> manager.set_api_key("example", "  not-a-real-key  ")
            >>> manager.get_api_key("example")
            'not-a-real-key'
            >>> manager.remove_api_key("example")
            True
        """
        config = self.load_config()
        if 'api_keys' not in config:
            config['api_keys'] = {}

        config['api_keys'][service] = api_key.strip()
        self.save_config(config)
        logging.info(f"API key saved for service: {service}")

    def remove_api_key(self, service: str) -> bool:
        """
        Remove API key for a specific service.

        Args:
            service: Service name.

        Returns:
            True if a key was removed, False if none was stored.

        Examples:
            >>> ConfigManager().remove_api_key("no-such-service")
            False
        """
        config = self.load_config()
        api_keys = config.get('api_keys', {})

        if service in api_keys:
            del api_keys[service]
            config['api_keys'] = api_keys
            self.save_config(config)
            logging.info(f"API key removed for service: {service}")
            return True
        return False

    def list_configured_services(self) -> list:
        """
        List all services with configured API keys.

        Returns:
            (list): Service names, in the order they were first stored.

        Examples:
            >>> manager = ConfigManager()
            >>> manager.set_api_key("example", "not-a-real-key")
            >>> "example" in manager.list_configured_services()
            True
            >>> manager.remove_api_key("example")
            True
        """
        config = self.load_config()
        api_keys = config.get('api_keys', {})
        return list(api_keys.keys())

    def get_config_info(self) -> Dict[str, Any]:
        """
        Get information about the configuration.

        Returns:
            (dict): ``config_directory``, ``config_file``, ``config_exists`` and
            ``configured_services``. The keys themselves are not included.

        Examples:
            >>> sorted(ConfigManager().get_config_info())
            ['config_directory', 'config_exists', 'config_file', 'configured_services']
        """
        return {
            'config_directory': str(self.config_dir),
            'config_file': str(self.config_file),
            'config_exists': self.config_file.exists(),
            'configured_services': self.list_configured_services()
        }


# Global configuration manager instance
_config_manager = None

def get_config_manager() -> ConfigManager:
    """
    Get the global configuration manager instance.

    Built on first call; the same object every time after.

    Returns:
        (ConfigManager): The shared instance.

    Examples:
        >>> get_config_manager() is get_config_manager()
        True
    """
    global _config_manager
    if _config_manager is None:
        _config_manager = ConfigManager()
    return _config_manager


def set_cas_api_key(api_key: str) -> Path:
    """
    Set CAS Common Chemistry API key for persistent storage

    [`CASCommonChem`][provesid.cascommonchem.CASCommonChem] uses the stored
    key from then on. The file's path is logged at INFO.

    Args:
        api_key: Your CAS API key

    Returns:
        (Path): The config file the key was written to.

    Note:
        This replaces any key already stored, without asking.

    Examples:
        >>> from provesid.config import set_cas_api_key
        >>> set_cas_api_key("your-cas-api-key-here")    # doctest: +SKIP
        PosixPath('/home/me/.config/provesid/config.json')
    """
    config_mgr = get_config_manager()
    config_mgr.set_api_key('cas', api_key)
    logger.info("CAS API key saved to %s; CASCommonChem() will use it", config_mgr.config_file)
    return config_mgr.config_file


def get_cas_api_key() -> Optional[str]:
    """
    Get the stored CAS API key.

    Only the config file is read; the ``CAS_API_KEY`` environment variable,
    which [`CASCommonChem`][provesid.cascommonchem.CASCommonChem] consults after this
    file, is not.

    Returns:
        The key, or None when none is stored.

    Examples:
        >>> key = get_cas_api_key()
        >>> key is None or isinstance(key, str)
        True
    """
    return get_config_manager().get_api_key('cas')


def remove_cas_api_key() -> bool:
    """
    Remove the stored CAS API key.

    Returns:
        True if a key was stored and is now gone, False if none was stored.

    Examples:
        >>> remove_cas_api_key()                          # doctest: +SKIP
        True
    """
    removed = get_config_manager().remove_api_key('cas')
    if removed:
        logger.info("CAS API key removed")
    else:
        logger.info("No CAS API key was configured")
    return removed


def show_config() -> Dict[str, Any]:
    """
    The configuration's location and which services have keys.

    The same as ``get_config_manager().get_config_info()``. The keys
    themselves are not included.

    Returns:
        (dict): ``config_directory``, ``config_file``, ``config_exists`` and
        ``configured_services``.

    Examples:
        >>> show_config()                                 # doctest: +SKIP
        {'config_directory': '/home/me/.config/provesid',
         'config_file': '/home/me/.config/provesid/config.json',
         'config_exists': True,
         'configured_services': ['cas']}
    """
    return get_config_manager().get_config_info()
