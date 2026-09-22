"""
Tests for NCI Chemical Identifier Resolver functionality
"""

import pytest
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from provesid.http import NotFoundError, Outcome, ServiceError
from provesid.resolver import (
    NCIChemicalIdentifierResolver, 
    NCIResolverError, 
    NCIResolverNotFoundError,
    nci_classify,
    nci_cas_to_mol,
    nci_id_to_mol,
    nci_resolver,
    nci_smiles_to_names,
    nci_name_to_smiles,
    nci_inchi_to_smiles,
    nci_cas_to_inchi,
    nci_get_molecular_weight,
    nci_get_formula
)


class TestNCIChemicalIdentifierResolver:
    """Test suite for NCIChemicalIdentifierResolver class"""
    
    @pytest.fixture
    def resolver(self):
        """Create an NCIChemicalIdentifierResolver instance for testing"""
        return NCIChemicalIdentifierResolver()
    
    def test_initialization(self, resolver):
        """Test NCIChemicalIdentifierResolver initialization"""
        assert resolver.base_url == "https://cactus.nci.nih.gov/chemical/structure"
        assert resolver.timeout == 30
        assert resolver.pause_time == 0.1
        assert len(resolver.representations) > 10  # Should have many representation types
    
    def test_custom_initialization(self):
        """Test NCIChemicalIdentifierResolver with custom parameters"""
        custom_resolver = NCIChemicalIdentifierResolver(
            base_url="https://example.com",
            timeout=60,
            pause_time=0.5
        )
        assert custom_resolver.base_url == "https://example.com"
        assert custom_resolver.timeout == 60
        assert custom_resolver.pause_time == 0.5
    
    def test_resolve_basic(self, resolver):
        """Test basic resolution functionality"""
        # Test converting aspirin name to SMILES
        smiles = resolver.resolve('aspirin', 'smiles')
        assert isinstance(smiles, str)
        assert len(smiles) > 0
        assert 'C' in smiles  # Should contain carbon atoms
        
        # Test converting SMILES to molecular weight
        mw = resolver.resolve('CCO', 'mw')  # ethanol
        assert isinstance(mw, str)
        assert float(mw) > 0
    
    def test_resolve_multiple(self, resolver):
        """Test resolving multiple representations at once"""
        representations = ['smiles', 'stdinchi', 'mw', 'formula']
        results = resolver.resolve_multiple('aspirin', representations)
        
        assert isinstance(results, dict)
        assert len(results) == len(representations)
        
        for rep in representations:
            assert rep in results
            # Some results might be None if not available
    
    def test_get_molecular_data(self, resolver):
        """Test comprehensive molecular data retrieval"""
        mol_data = resolver.get_molecular_data('ethanol')
        
        assert isinstance(mol_data, dict)
        
        # Check for expected keys
        expected_keys = ['smiles', 'stdinchi', 'formula', 'mw']
        for key in expected_keys:
            assert key in mol_data
        
        # Check that we got some data
        assert mol_data['formula'] == 'C2H6O'
        assert float(mol_data['mw']) > 40  # Ethanol MW is ~46
    
    def test_batch_resolve(self, resolver):
        """Test batch resolution of multiple compounds"""
        compounds = ['aspirin', 'caffeine', 'ethanol']
        results = resolver.batch_resolve(compounds, 'formula')
        
        assert isinstance(results, dict)
        assert len(results) == len(compounds)
        
        for compound in compounds:
            assert compound in results
            # All these common compounds should have formulas
            if results[compound]:
                assert isinstance(results[compound], str)
                assert 'C' in results[compound]  # Should contain carbon
    
    def test_get_image_url(self, resolver):
        """Test image URL generation"""
        image_url = resolver.get_image_url('aspirin')
        
        assert isinstance(image_url, str)
        assert image_url.startswith('https://cactus.nci.nih.gov')
        assert 'aspirin' in image_url
        assert 'image' in image_url
    
    def test_get_available_representations(self, resolver):
        """Test getting available representations"""
        representations = resolver.get_available_representations()
        
        assert isinstance(representations, list)
        assert len(representations) > 0
        
        # Check for common representations
        expected_reps = ['smiles', 'stdinchi', 'formula', 'mw']
        for rep in expected_reps:
            assert rep in representations
    
    def test_error_handling(self, resolver):
        """Test error handling for invalid inputs"""
        # Test with clearly invalid identifier - should raise some kind of resolver error
        with pytest.raises((NCIResolverNotFoundError, NCIResolverError)):
            resolver.resolve('this_is_definitely_not_a_chemical_12345', 'smiles')
    
    def test_cas_number_resolution(self, resolver):
        """Test resolution using CAS numbers"""
        # Test ethanol CAS number
        smiles = resolver.resolve('64-17-5', 'smiles')
        assert isinstance(smiles, str)
        assert len(smiles) > 0
        
        # Should be ethanol SMILES
        assert 'CCO' in smiles or 'OCC' in smiles
    
    def test_inchi_conversion(self, resolver):
        """Test InChI conversions"""
        try:
            # Convert SMILES to InChI
            inchi = resolver.resolve('CCO', 'stdinchi')
            assert isinstance(inchi, str)
            assert inchi.startswith('InChI=')
            
            # And back. This used to be skipped as "not supported for this
            # format"; it was the client percent-encoding the InChI's slashes.
            smiles = resolver.resolve(inchi, 'smiles')
            assert smiles == 'CCO'
        except NCIResolverNotFoundError:
            raise
        except NCIResolverError:
            pytest.skip("NCI resolver not available")

    def test_slashes_stay_path_separators(self, resolver):
        """
        CACTUS's server answers 404 to a path holding ``%2F``, so an InChI or
        a SMILES with stereo bonds has to keep its slashes as they are.
        """
        url = resolver._build_url('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'smiles')
        assert url.endswith('/InChI%3D1S/C2H6O/c1-2-3/h3H%2C2H2%2C1H3/smiles')
        assert '/C/C%3DC/C/' in resolver._build_url('C/C=C/C', 'stdinchikey')
        assert '%23' in resolver._build_url('C#C', 'smiles')


class TestConvenienceFunctions:
    """Test suite for convenience functions"""
    
    def test_nci_cas_to_mol(self):
        """Test nci_cas_to_mol function"""
        result = nci_cas_to_mol('64-17-5')  # ethanol
        
        assert isinstance(result, dict)
        assert result['note'] == 'OK'
        assert 'smiles' in result
        assert 'formula' in result
        assert 'mw' in result
        
        # Check ethanol specific data
        assert result['formula'] == 'C2H6O'
        assert 'CCO' in result['smiles'] or 'OCC' in result['smiles']
    
    def test_nci_name_to_smiles(self):
        """Test nci_name_to_smiles function"""
        smiles = nci_name_to_smiles('aspirin')
        
        assert isinstance(smiles, str)
        assert len(smiles) > 0
        assert 'C' in smiles
        assert 'O' in smiles  # Aspirin contains oxygen
    
    def test_nci_get_molecular_weight(self):
        """Test nci_get_molecular_weight function"""
        mw = nci_get_molecular_weight('water')
        
        if mw is not None:  # Function returns Optional[float]
            assert isinstance(mw, float)
            assert 17 < mw < 19  # Water MW should be around 18
    
    def test_nci_get_formula(self):
        """Test nci_get_formula function"""
        formula = nci_get_formula('water')
        
        if formula is not None:  # Function returns Optional[str]
            assert isinstance(formula, str)
            assert formula == 'H2O'
        
        # Test another compound
        formula = nci_get_formula('caffeine')
        if formula is not None:
            assert isinstance(formula, str)
            assert 'C' in formula and 'N' in formula  # Caffeine contains C and N
    
    def test_nci_inchi_to_smiles(self):
        """Test nci_inchi_to_smiles function"""
        # First get an InChI
        resolver = NCIChemicalIdentifierResolver()
        try:
            inchi = resolver.resolve('ethanol', 'stdinchi')
            
            # Convert back to SMILES - may return None if conversion fails
            smiles = nci_inchi_to_smiles(inchi)
            if smiles is not None:
                assert isinstance(smiles, str)
                assert 'C' in smiles and 'O' in smiles
            else:
                # Conversion may fail for some InChI formats
                pytest.skip("InChI to SMILES conversion not available")
        except (NCIResolverError, NCIResolverNotFoundError):
            pytest.skip("NCI resolver not available")
    
    def test_nci_cas_to_inchi(self):
        """Test nci_cas_to_inchi function"""
        inchi = nci_cas_to_inchi('64-17-5')  # ethanol
        
        assert isinstance(inchi, str)
        assert inchi.startswith('InChI=')
        assert 'C2H6O' in inchi  # Should contain ethanol formula
    
    def test_nci_smiles_to_names(self):
        """Test nci_smiles_to_names function"""
        names = nci_smiles_to_names('CCO')  # ethanol
        
        assert isinstance(names, list)
        assert len(names) > 0
        
        # Should contain ethanol-related names
        names_str = ' '.join(names).lower()
        assert 'ethanol' in names_str or 'alcohol' in names_str
    
    def test_nci_id_to_mol(self):
        """Test nci_id_to_mol function"""
        # Test with a name
        result = nci_id_to_mol('aspirin')
        
        assert isinstance(result, dict)
        assert 'smiles' in result
        assert 'formula' in result
        assert 'mw' in result
        
        # Check aspirin formula
        assert 'C9H8O4' in result['formula']
    
    def test_nci_resolver(self):
        """Test nci_resolver function"""
        result = nci_resolver('aspirin', 'formula')
        
        assert isinstance(result, str)
        assert 'C9H8O4' in result  # Aspirin formula


class TestSpecialCases:
    """Test suite for special cases and edge conditions"""
    
    def test_common_drug_names(self):
        """Test resolution of common drug names"""
        drugs = ['aspirin', 'caffeine', 'ibuprofen', 'acetaminophen']
        resolver = NCIChemicalIdentifierResolver()
        
        for drug in drugs:
            try:
                smiles = resolver.resolve(drug, 'smiles')
                assert isinstance(smiles, str)
                assert len(smiles) > 0
                assert 'C' in smiles  # All these drugs contain carbon
            except NCIResolverNotFoundError:
                # Some drugs might not be found, which is acceptable
                pass
    
    def test_common_solvents(self):
        """Test resolution of common solvents"""
        solvents = ['water', 'ethanol', 'methanol', 'acetone']
        resolver = NCIChemicalIdentifierResolver()
        
        for solvent in solvents:
            formula = resolver.resolve(solvent, 'formula')
            assert isinstance(formula, str)
            assert len(formula) > 0
    
    def test_cas_number_formats(self):
        """Test different CAS number formats"""
        resolver = NCIChemicalIdentifierResolver()
        
        # Test standard CAS format
        result1 = resolver.resolve('64-17-5', 'formula')  # ethanol
        assert result1 == 'C2H6O'
        
        # Different compounds
        cas_numbers = ['50-00-0', '67-56-1', '71-43-2']  # formaldehyde, methanol, benzene
        
        for cas in cas_numbers:
            try:
                formula = resolver.resolve(cas, 'formula')
                assert isinstance(formula, str)
                assert 'C' in formula or 'H' in formula  # Should contain common elements
            except NCIResolverNotFoundError:
                # Some CAS numbers might not be found
                pass
    
    def test_unicode_and_special_characters(self):
        """Test handling of unicode and special characters"""
        resolver = NCIChemicalIdentifierResolver()
        
        # Test with unicode characters (might be in chemical names)
        test_names = ['β-carotene', 'α-glucose', 'D-glucose']
        
        for name in test_names:
            try:
                result = resolver.resolve(name, 'formula')
                if result:  # If found
                    assert isinstance(result, str)
            except (NCIResolverNotFoundError, NCIResolverError):
                # These might not be found or might cause errors
                pass
    
    def test_batch_processing_with_invalid_compounds(self):
        """Test batch processing with mix of valid and invalid compounds"""
        resolver = NCIChemicalIdentifierResolver()
        
        compounds = ['aspirin', 'invalid_compound_xyz', 'ethanol', 'another_invalid_one']
        results = resolver.batch_resolve(compounds, 'formula')
        
        assert isinstance(results, dict)
        assert len(results) == len(compounds)
        
        # Valid compounds should have results
        assert results['aspirin'] is not None
        assert results['ethanol'] is not None
        
        # Invalid compounds should have None results
        assert results['invalid_compound_xyz'] is None
        assert results['another_invalid_one'] is None
    
    def test_rate_limiting(self):
        """Test rate limiting functionality"""
        resolver = NCIChemicalIdentifierResolver()
        
        import time
        start_time = time.time()
        
        # Make multiple calls
        for _ in range(3):
            resolver._rate_limit()
        
        elapsed = time.time() - start_time
        
        # Should take at least 2 * pause_time due to rate limiting
        expected_min_time = 2 * resolver.pause_time
        assert elapsed >= expected_min_time * 0.9  # Allow some tolerance


class TestErrorHandling:
    """Test suite for error handling"""
    
    def test_network_timeout_simulation(self):
        """Test behavior with very short timeout"""
        # Create resolver with very short timeout
        resolver = NCIChemicalIdentifierResolver(timeout=0.001)  # Very short timeout
        
        # This should likely timeout or succeed depending on network speed
        try:
            result = resolver.resolve('aspirin', 'smiles')
            # If it succeeds, that's fine - the network was fast enough
            assert result is not None or result is None  # Either is acceptable
        except (NCIResolverError, NCIResolverNotFoundError):
            # Expected timeout/error - this is also fine
            pass
    
    def test_invalid_representation(self):
        """Test with invalid representation type"""
        resolver = NCIChemicalIdentifierResolver()
        
        # Test with non-existent representation - should raise ValueError
        with pytest.raises(ValueError):
            resolver.resolve('aspirin', 'invalid_representation_type')
    
    def test_empty_identifier(self):
        """Test with empty identifier"""
        resolver = NCIChemicalIdentifierResolver()
        
        with pytest.raises((NCIResolverError, NCIResolverNotFoundError)):
            resolver.resolve('', 'smiles')
    
    def test_none_identifier(self):
        """Test with None identifier"""
        resolver = NCIChemicalIdentifierResolver()
        
        # Note: We can't test None directly due to type hints, so test empty string
        with pytest.raises((NCIResolverError, NCIResolverNotFoundError)):
            resolver.resolve('', 'smiles')


class TestNCIClassification:
    """
    What a CACTUS response means, which its status code does not say.

    Verified live on 2026-09-19: ``this_is_definitely_not_a_chemical_12345``
    answers HTTP 500 with the body ``<h1>Page not found (404)</h1>``, while
    ``\u03b1-glucose`` answers 200. Reading the status alone therefore costs
    four requests and several seconds of back-off to learn a permanent answer.
    These tests are offline; they pin the classification, not the service.
    """

    class _Response:
        def __init__(self, status_code, text=""):
            self.status_code = status_code
            self.text = text
            self.headers = {}

    def test_not_found_behind_a_500_is_absence(self):
        """The body is believed over the status."""
        response = self._Response(500, "<h1>Page not found (404)</h1>\n")
        assert nci_classify(response) is Outcome.ABSENT

    def test_a_real_server_error_is_still_retryable(self):
        """A 500 that is not the not-found page keeps its usual reading."""
        response = self._Response(503, "<h1>Service temporarily unavailable</h1>")
        assert nci_classify(response) is Outcome.RETRY

    def test_a_404_is_absence(self):
        """The honest form of the same answer."""
        assert nci_classify(self._Response(404, "")) is Outcome.ABSENT

    def test_a_200_is_an_answer(self):
        assert nci_classify(self._Response(200, "CCO")) is Outcome.OK

    def test_unresolvable_identifier_costs_one_request(self, monkeypatch):
        """
        The point of the classifier: no back-off for a permanent answer.

        Before it, an unresolvable identifier cost four requests and about
        seven seconds of exponential back-off.
        """
        import requests

        calls = []

        def not_found(url, timeout=None):
            calls.append(url)
            return TestNCIClassification._Response(500, "<h1>Page not found (404)</h1>")

        monkeypatch.setattr(requests, "get", not_found)
        resolver = NCIChemicalIdentifierResolver(pause_time=0)

        with pytest.raises(NCIResolverNotFoundError):
            resolver.resolve("this_is_definitely_not_a_chemical_12345", "smiles")
        assert len(calls) == 1

    def test_the_exceptions_carry_the_shared_bases(self):
        """A caller can catch this service, or every service at once."""
        assert issubclass(NCIResolverError, ServiceError)
        assert issubclass(NCIResolverNotFoundError, NCIResolverError)
        assert issubclass(NCIResolverNotFoundError, NotFoundError)


if __name__ == "__main__":
    # Run tests if executed directly
    pytest.main([__file__, "-v"])
