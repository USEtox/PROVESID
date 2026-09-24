"""
provesid.config's module functions return what they used to print, and log
it; they never write to stdout (dev-principle 8).
"""

import logging

import pytest

from provesid import config


@pytest.fixture
def manager(tmp_path, monkeypatch):
    """A config manager in a temporary directory, in place of the shared one."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.setenv("APPDATA", str(tmp_path))
    manager = config.ConfigManager()
    monkeypatch.setattr(config, "_config_manager", manager)
    return manager


@pytest.mark.unit
def test_set_returns_the_file_and_logs(manager, capsys, caplog):
    with caplog.at_level(logging.INFO, logger="provesid.config"):
        path = config.set_cas_api_key("not-a-real-key")

    assert path == manager.config_file
    assert config.get_cas_api_key() == "not-a-real-key"
    assert str(path) in caplog.text
    assert capsys.readouterr().out == ""


@pytest.mark.unit
def test_remove_says_whether_there_was_a_key(manager, capsys):
    config.set_cas_api_key("not-a-real-key")

    assert config.remove_cas_api_key() is True
    assert config.remove_cas_api_key() is False
    assert config.get_cas_api_key() is None
    assert capsys.readouterr().out == ""


@pytest.mark.unit
def test_show_config_returns_the_info(manager, capsys):
    config.set_cas_api_key("not-a-real-key")

    info = config.show_config()

    assert info["config_file"] == str(manager.config_file)
    assert info["configured_services"] == ["cas"]
    assert "not-a-real-key" not in str(info)
    assert capsys.readouterr().out == ""
