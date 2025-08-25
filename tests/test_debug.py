"""Simple test to check pytest collection"""

import pytest
from provesid.pubchem import PubChemAPI

class TestPubChemSimple:
    
    def test_basic(self):
        """Basic test"""
        api = PubChemAPI()
        assert api.base_url == "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        
    def test_another(self):
        """Another test"""
        assert True
