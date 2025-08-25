"""
Pytest configuration and shared fixtures for PROVESID tests.
"""

import pytest
import requests
import time
from unittest.mock import Mock, patch
import os
import sys

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


# Pytest markers for test categorization
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "unit: mark test as a unit test"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test"
    )
    config.addinivalue_line(
        "markers", "api: mark test as requiring API access"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "pubchem: mark test as PubChem specific"
    )
    config.addinivalue_line(
        "markers", "cas: mark test as CAS Common Chemistry specific"
    )
    config.addinivalue_line(
        "markers", "chebi: mark test as ChEBI specific"
    )
    config.addinivalue_line(
        "markers", "classyfire: mark test as ClassyFire specific"
    )
    config.addinivalue_line(
        "markers", "opsin: mark test as OPSIN specific"
    )
    config.addinivalue_line(
        "markers", "nci: mark test as NCI Resolver specific"
    )


@pytest.fixture(scope="session")
def api_timeout():
    """Default timeout for API calls in tests."""
    return 30


@pytest.fixture(scope="session")
def test_compounds():
    """Common test compounds used across multiple test modules."""
    return {
        'water': {
            'name': 'water',
            'cas': '7732-18-5',
            'smiles': 'O',
            'cid': 962,
            'formula': 'H2O',
            'chebi_id': 'CHEBI:15377'
        },
        'ethanol': {
            'name': 'ethanol',
            'cas': '64-17-5',
            'smiles': 'CCO',
            'cid': 702,
            'formula': 'C2H6O',
            'chebi_id': 'CHEBI:16236'
        },
        'aspirin': {
            'name': 'aspirin',
            'cas': '50-78-2',
            'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
            'cid': 2244,
            'formula': 'C9H8O4',
            'chebi_id': 'CHEBI:15365'
        },
        'caffeine': {
            'name': 'caffeine',
            'cas': '58-08-2',
            'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'cid': 2519,
            'formula': 'C8H10N4O2',
            'chebi_id': 'CHEBI:27732'
        }
    }


@pytest.fixture
def mock_response():
    """Create a mock HTTP response for testing."""
    def _mock_response(status_code=200, json_data=None, text=None, headers=None):
        mock = Mock()
        mock.status_code = status_code
        mock.headers = headers or {}
        
        if json_data is not None:
            mock.json.return_value = json_data
        else:
            mock.json.side_effect = ValueError("No JSON object could be decoded")
        
        mock.text = text or ""
        mock.content = (text or "").encode('utf-8')
        return mock
    
    return _mock_response


@pytest.fixture
def rate_limiter():
    """Rate limiter to prevent overwhelming APIs during tests."""
    last_call = {'time': 0}
    min_interval = 0.1  # 100ms between calls
    
    def _rate_limit():
        current_time = time.time()
        time_since_last = current_time - last_call['time']
        if time_since_last < min_interval:
            time.sleep(min_interval - time_since_last)
        last_call['time'] = time.time()
    
    return _rate_limit


@pytest.fixture
def skip_if_no_internet():
    """Skip test if no internet connection is available."""
    try:
        requests.get('http://httpbin.org/status/200', timeout=5)
        return False
    except requests.RequestException:
        pytest.skip("No internet connection available")


@pytest.fixture
def skip_if_api_unavailable():
    """Skip test if API endpoints are unavailable."""
    def _check_api(url, timeout=5):
        try:
            response = requests.get(url, timeout=timeout)
            if response.status_code >= 500:
                pytest.skip(f"API at {url} is unavailable (status: {response.status_code})")
            return True
        except requests.RequestException:
            pytest.skip(f"API at {url} is not reachable")
    
    return _check_api


class APITestHelper:
    """Helper class for API testing utilities."""
    
    @staticmethod
    def is_valid_json_response(response):
        """Check if response is valid JSON."""
        try:
            response.json()
            return True
        except (ValueError, AttributeError):
            return False
    
    @staticmethod
    def assert_valid_compound_data(data, compound_info=None):
        """Assert that compound data has expected structure."""
        assert isinstance(data, dict)
        
        if compound_info:
            # Check specific compound properties if provided
            if 'formula' in compound_info and 'molecularFormula' in data:
                assert compound_info['formula'] in data['molecularFormula']
    
    @staticmethod
    def mock_api_error(error_type='timeout'):
        """Create mock API error for testing."""
        if error_type == 'timeout':
            return requests.Timeout("Request timed out")
        elif error_type == 'connection':
            return requests.ConnectionError("Connection failed")
        elif error_type == 'http':
            return requests.HTTPError("HTTP error")
        else:
            return requests.RequestException("Generic request error")


@pytest.fixture
def api_helper():
    """Provide API testing helper utilities."""
    return APITestHelper()


# Test data fixtures for specific APIs
@pytest.fixture
def pubchem_test_data():
    """PubChem-specific test data."""
    return {
        'valid_cids': [2244, 702, 962],  # aspirin, ethanol, water
        'invalid_cids': [999999999, -1, 0],
        'valid_names': ['aspirin', 'ethanol', 'water'],
        'invalid_names': ['xyz123notachemical', '', '!@#$%'],
        'valid_smiles': ['CCO', 'O', 'CC(=O)OC1=CC=CC=C1C(=O)O'],
        'invalid_smiles': ['invalid_smiles', 'XYZ123', '']
    }


@pytest.fixture
def cas_test_data():
    """CAS Common Chemistry test data."""
    return {
        'valid_cas': ['7732-18-5', '64-17-5', '50-78-2'],
        'invalid_cas': ['invalid-cas', '123-45-6789', ''],
        'valid_names': ['water', 'ethanol', 'aspirin'],
        'invalid_names': ['xyz123notachemical', '', '!@#$%']
    }


@pytest.fixture
def chebi_test_data():
    """ChEBI test data."""
    return {
        'valid_ids': ['CHEBI:15377', 'CHEBI:16236', 'CHEBI:15365'],
        'invalid_ids': ['CHEBI:999999', 'INVALID:123', ''],
        'valid_names': ['water', 'ethanol', 'aspirin'],
        'invalid_names': ['xyz123notachemical', '', '!@#$%']
    }


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers based on test file names."""
    for item in items:
        # Add markers based on test file names
        if "test_pubchem" in item.nodeid:
            item.add_marker(pytest.mark.pubchem)
        elif "test_cas" in item.nodeid:
            item.add_marker(pytest.mark.cas)
        elif "test_chebi" in item.nodeid:
            item.add_marker(pytest.mark.chebi)
        elif "test_classyfire" in item.nodeid:
            item.add_marker(pytest.mark.classyfire)
        elif "test_opsin" in item.nodeid:
            item.add_marker(pytest.mark.opsin)
        elif "test_nci" in item.nodeid:
            item.add_marker(pytest.mark.nci)
        
        # Mark tests that make API calls
        if any(marker in item.keywords for marker in ['integration', 'api']):
            item.add_marker(pytest.mark.api)
