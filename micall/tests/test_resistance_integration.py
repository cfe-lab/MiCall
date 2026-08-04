"""
Integration tests for resistance reporting pipeline.

These tests exercise the complete flow from amino acid counts through resistance
calculation to PDF report generation, catching configuration and integration issues
that unit tests miss.
"""
import pytest
from io import StringIO, BytesIO
from csv import DictReader

from micall.resistance.resistance import (
    write_resistance,
    load_asi,
    read_aminos,
    filter_aminos,
    REPORTED_REGIONS
)
from micall.resistance.genreport import gen_report, read_config


@pytest.fixture
def asi_algorithms():
    """Load ASI algorithms for testing."""
    return load_asi()


def create_hiv_amino_csv(regions_data):
    """Helper to create amino CSV for HIV regions.
    
    Args:
        regions_data: dict of {region_name: [(position, amino_dict, coverage)]}
    
    Returns:
        StringIO with formatted CSV
    """
    lines = [
        'seed,region,q-cutoff,query.nuc.pos,refseq.aa.pos,'
        'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,'
        'partial,del,ins,clip,g2p_overlap,coverage'
    ]
    
    for region, positions in regions_data.items():
        for pos, amino_counts, coverage in positions:
            amino_values = [amino_counts.get(aa, 0) for aa in 'ACDEFGHIKLMNPQRSTVWY*X']
            amino_str = ','.join(str(v) for v in amino_values)
            lines.append(
                f'HIV1-B-FR-K03455-seed,{region},15,{pos},{pos},'
                f'{amino_str},0,0,0,0,0,{coverage}'
            )
    
    return StringIO('\n'.join(lines))


class TestHIVCARegionIntegration:
    """Test the complete pipeline for HIV CA region."""
    
    def test_ca_region_through_resistance_pipeline(self, asi_algorithms):
        """Test CA region flows through write_resistance without errors.
        
        This test ensures:
        1. CA amino data is read correctly
        2. CA is processed by resistance algorithm
        3. CA resistance calls are written to CSV
        4. No KeyError or other exceptions occur
        """
        # Create amino CSV with CA data covering all key positions
        # CA key positions: [56, 57, 66, 67, 70, 74, 105, 107]
        # Need to cover up to position 107
        ca_data = [
            *[(i, {'A': 100}, 100) for i in range(1, 108)],
        ]
        
        amino_csv = create_hiv_amino_csv({'CA': ca_data})
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        
        # Read aminos (need DictReader)
        amino_csv.seek(0)
        aminos = list(read_aminos(DictReader(amino_csv), min_fraction=0.05, min_coverage=100,
                                  algorithms=asi_algorithms))
        
        # Filter aminos (adds missing regions)
        filtered_aminos = filter_aminos(aminos, asi_algorithms)
        
        # Write resistance
        write_resistance(filtered_aminos,
                        resistance_csv,
                        mutations_csv,
                        algorithms=asi_algorithms)
        
        # Verify CA appears in resistance output
        resistance_csv.seek(0)
        resistance_lines = list(DictReader(resistance_csv))
        
        # Should have CA entries
        ca_entries = [r for r in resistance_lines if r['region'] == 'CA']
        assert len(ca_entries) > 0, "CA region should appear in resistance output"
        
        # Verify drug codes are present (Lenacapavir for CA)
        ca_drugs = {r['drug'] for r in ca_entries}
        assert 'LEN' in ca_drugs, "Lenacapavir (LEN) should be in CA resistance calls"
    
    def test_ca_region_through_complete_pipeline(self, asi_algorithms):
        """Test CA region through full pipeline including PDF generation.
        
        This is the ultimate integration test - if this passes, CA is fully integrated.
        """
        # Create comprehensive amino CSV with CA and other regions
        # Must cover key positions: CA max=107, PR max=90, RT max=348
        regions_data = {
            'CA': [(i, {'A': 100}, 100) for i in range(1, 108)],
            'PR': [(i, {'K': 100}, 100) for i in range(1, 91)],
            'RT': [(i, {'C': 100}, 100) for i in range(1, 349)],
        }
        
        amino_csv = create_hiv_amino_csv(regions_data)
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        
        # Step 1: Read and filter aminos
        amino_csv.seek(0)
        aminos = list(read_aminos(DictReader(amino_csv), min_fraction=0.05, min_coverage=100,
                                  algorithms=asi_algorithms))
        filtered_aminos = filter_aminos(aminos, asi_algorithms)
        
        # Step 2: Generate resistance calls
        write_resistance(filtered_aminos,
                        resistance_csv,
                        mutations_csv,
                        algorithms=asi_algorithms)
        
        # Step 3: Generate PDF report (using BytesIO to avoid file system)
        resistance_csv.seek(0)
        mutations_csv.seek(0)
        pdf_output = BytesIO()
        
        # This should not raise ValueError: Unknown drug codes or KeyError: 'CA'
        try:
            gen_report(resistance_csv, mutations_csv, pdf_output,
                      sample_name='TEST_CA_SAMPLE', git_version='test')
        except (KeyError, ValueError) as e:
            pytest.fail(f"PDF generation failed with {type(e).__name__}: {e}. "
                       f"This indicates configuration mismatch between resistance.py "
                       f"and genreport.yaml")
        
        # Verify PDF was generated (has content)
        assert pdf_output.tell() > 0, "PDF should have been written"


class TestHIVAllRegionsIntegration:
    """Test complete pipeline with all HIV regions."""
    
    def test_all_hiv_regions_end_to_end(self, asi_algorithms):
        """Test all HIV regions (PR, RT, IN, CA) through complete pipeline.
        
        This test ensures no region causes pipeline failures and all are
        properly integrated with genreport.yaml configuration.
        """
        # Create data for all HIV regions covering their key positions
        # Key position ranges: PR max=90, RT max=348, IN max=263, CA max=107
        regions_data = {
            'PR': [(i, {'P': 100}, 100) for i in range(1, 91)],
            'RT': [(i, {'M': 100}, 100) for i in range(1, 349)],
            'INT': [(i, {'E': 100}, 100) for i in range(1, 264)],
            'CA': [(i, {'A': 100}, 100) for i in range(1, 108)],
        }
        
        amino_csv = create_hiv_amino_csv(regions_data)
        resistance_csv = StringIO()
        mutations_csv = StringIO()
        
        # Full pipeline
        amino_csv.seek(0)
        aminos = list(read_aminos(DictReader(amino_csv), min_fraction=0.05, min_coverage=100,
                                  algorithms=asi_algorithms))
        filtered_aminos = filter_aminos(aminos, asi_algorithms)
        write_resistance(filtered_aminos, resistance_csv, mutations_csv,
                        algorithms=asi_algorithms)
        
        # Verify all regions appear in output
        resistance_csv.seek(0)
        resistance_lines = list(DictReader(resistance_csv))
        
        output_regions = {r['coord_region'] for r in resistance_lines}
        expected_regions = {'PR', 'RT', 'INT', 'CA'}
        
        missing = expected_regions - output_regions
        assert not missing, f"Missing regions in output: {missing}"
        
        # Test PDF generation
        resistance_csv.seek(0)
        mutations_csv.seek(0)
        pdf_output = BytesIO()
        
        gen_report(resistance_csv, mutations_csv, pdf_output,
                  sample_name='TEST_ALL_REGIONS', git_version='test')
        
        assert pdf_output.tell() > 0, "PDF should have been generated"


class TestConfigurationConsistency:
    """Test configuration consistency between modules."""
    
    def test_reported_regions_have_algorithm_support(self, asi_algorithms):
        """Verify regions in REPORTED_REGIONS have algorithm support.
        
        This test ensures we're not trying to report resistance for regions
        that don't have scoring rules in the algorithm.
        """
        hiv_algorithm = asi_algorithms[None]
        algorithm_regions = set(hiv_algorithm.gene_def.keys())
        
        # Map INT -> IN for comparison
        if 'INT' in REPORTED_REGIONS and 'IN' in algorithm_regions:
            algorithm_regions.add('INT')
        
        # Filter to HIV regions only
        hiv_reported = {r for r in REPORTED_REGIONS 
                       if r not in ('NS3', 'NS5a', 'NS5b')}
        
        unsupported = hiv_reported - algorithm_regions
        assert not unsupported, (
            f"HIV regions in REPORTED_REGIONS lack algorithm support: {unsupported}"
        )
    
    def test_algorithm_regions_processed_correctly(self, asi_algorithms):
        """Test that all algorithm regions can be processed.
        
        Ensures get_algorithm_regions and filter_aminos handle all regions
        from the algorithm without errors.
        """
        from micall.resistance.resistance import get_algorithm_regions, create_empty_aminos
        
        hiv_algorithm = asi_algorithms[None]
        regions = get_algorithm_regions(hiv_algorithm)
        
        # Should be able to create empty aminos for each region
        for region in regions:
            try:
                empty = create_empty_aminos(region, None, 'HIV1B-seed', asi_algorithms)
                assert empty.region == region
                assert len(empty.aminos) > 0
            except Exception as e:
                pytest.fail(f"Failed to create empty aminos for region {region}: {e}")


class TestRegressionPrevention:
    """Tests that prevent specific past bugs from recurring."""
    
    def test_ca_in_reported_regions_and_genreport(self):
        """Regression test: Ensure CA is in both REPORTED_REGIONS and genreport.yaml.
        
        This was the bug that caused the microtest failure. This test ensures
        it can never happen again.
        """
        # CA must be in REPORTED_REGIONS
        assert 'CA' in REPORTED_REGIONS, (
            "CA must be in REPORTED_REGIONS to enable CA resistance reporting"
        )
        
        # CA must be in genreport.yaml
        config = read_config('test')
        report_template = config[0]
        known_regions = report_template.virus_config['known_regions']
        
        assert 'CA' in known_regions, (
            "CA must be in genreport.yaml known_regions to process CA resistance data"
        )
    
    def test_lenacapavir_in_genreport(self):
        """Regression test: Ensure Lenacapavir (LEN) is configured in genreport.yaml.
        
        This was the second bug - LEN was in algorithm but not in known_drugs.
        """
        config = read_config('test')
        report_template = config[0]
        known_drugs = report_template.virus_config['known_drugs']
        
        # Extract all drug codes
        all_drug_codes = {
            drug_code
            for drug_class in known_drugs.values()
            for drug_code, drug_name in drug_class
        }
        
        assert 'LEN' in all_drug_codes, (
            "Lenacapavir (LEN) must be in genreport.yaml known_drugs for CA resistance reporting"
        )
        
        # Verify it's in CAI drug class
        assert 'CAI' in known_drugs, "CAI (Capsid Inhibitor) drug class must exist"
        cai_drugs = {drug_code for drug_code, _ in known_drugs['CAI']}
        assert 'LEN' in cai_drugs, "Lenacapavir must be in CAI drug class"
