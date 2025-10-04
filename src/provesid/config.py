"""
Configuration management for PROVESID API keys and settings.
Provides persistent storage for API keys and user preferences.
"""

import os
import json
import logging
from pathlib import Path
from typing import Optional, Dict, Any


class ConfigManager:
    """Manages persistent configuration for PROVESID"""
    
    def __init__(self):
        """Initialize configuration manager with default paths"""
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
        """Load configuration from file"""
        if not self.config_file.exists():
            return {}
        
        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            logging.warning(f"Error loading config from {self.config_file}: {e}")
            return {}
    
    def save_config(self, config: Dict[str, Any]):
        """Save configuration to file"""
        try:
            with open(self.config_file, 'w', encoding='utf-8') as f:
                json.dump(config, f, indent=2, ensure_ascii=False)
        except Exception as e:
            logging.error(f"Error saving config to {self.config_file}: {e}")
    
    def get_api_key(self, service: str) -> Optional[str]:
        """Get API key for a specific service"""
        config = self.load_config()
        api_keys = config.get('api_keys', {})
        return api_keys.get(service)
    
    def set_api_key(self, service: str, api_key: str):
        """Set API key for a specific service"""
        config = self.load_config()
        if 'api_keys' not in config:
            config['api_keys'] = {}
        
        config['api_keys'][service] = api_key.strip()
        self.save_config(config)
        logging.info(f"API key saved for service: {service}")
    
    def remove_api_key(self, service: str) -> bool:
        """Remove API key for a specific service"""
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
        """List all services with configured API keys"""
        config = self.load_config()
        api_keys = config.get('api_keys', {})
        return list(api_keys.keys())
    
    def get_config_info(self) -> Dict[str, Any]:
        """Get information about the configuration"""
        return {
            'config_directory': str(self.config_dir),
            'config_file': str(self.config_file),
            'config_exists': self.config_file.exists(),
            'configured_services': self.list_configured_services()
        }


# Global configuration manager instance
_config_manager = None

def get_config_manager() -> ConfigManager:
    """Get the global configuration manager instance"""
    global _config_manager
    if _config_manager is None:
        _config_manager = ConfigManager()
    return _config_manager


def set_cas_api_key(api_key: str):
    """
    Set CAS Common Chemistry API key for persistent storage
    
    Args:
        api_key: Your CAS API key
        
    Example:
        >>> from provesid.config import set_cas_api_key
        >>> set_cas_api_key("your-cas-api-key-here")
        >>> # Now CASCommonChem() will automatically use this key
    """
    config_mgr = get_config_manager()
    config_mgr.set_api_key('cas', api_key)
    print(f"✅ CAS API key saved to: {config_mgr.config_file}")
    print("ℹ️  CASCommonChem() will now automatically use this key")


def get_cas_api_key() -> Optional[str]:
    """Get the stored CAS API key"""
    return get_config_manager().get_api_key('cas')


def remove_cas_api_key():
    """Remove the stored CAS API key"""
    config_mgr = get_config_manager()
    if config_mgr.remove_api_key('cas'):
        print("✅ CAS API key removed")
    else:
        print("ℹ️  No CAS API key was configured")


def show_config():
    """Show current configuration information"""
    config_mgr = get_config_manager()
    info = config_mgr.get_config_info()
    
    print("PROVESID Configuration:")
    print(f"  Config directory: {info['config_directory']}")
    print(f"  Config file: {info['config_file']}")
    print(f"  Config exists: {info['config_exists']}")
    print(f"  Configured services: {', '.join(info['configured_services']) if info['configured_services'] else 'None'}")