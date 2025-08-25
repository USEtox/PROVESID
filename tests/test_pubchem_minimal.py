"""
Minimal PubChem tests for debugging
"""

import pytest
from unittest.mock import patch, Mock

from provesid.pubchem import (
    PubChemAPI,
    PubChemError,
    PubChemNotFoundError,
    PubChemTimeoutError,
    PubChemServerError,
    CompoundProperties,
    OutputFormat,
    Domain,
    CompoundDomainNamespace,
    Operation
)


@pytest.mark.unit
class TestPubChemAPI:
    """Test suite for PubChemAPI class"""
    
    @pytest.fixture
    def api(self):
        """Create a PubChemAPI instance for testing"""
        return PubChemAPI()
    
    def test_initialization(self, api):
        """Test PubChemAPI initialization"""
        assert api.base_url == "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        assert api.pause_time == 0.2
    
    def test_custom_initialization(self):
        """Test PubChemAPI with custom parameters"""
        custom_api = PubChemAPI(
            base_url="https://example.com",
            pause_time=0.5
        )
        assert custom_api.base_url == "https://example.com"
        assert custom_api.pause_time == 0.5
