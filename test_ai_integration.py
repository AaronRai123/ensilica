"""
Test script for the AI integration in Ensilica
"""
import os
import json
import logging
from modules.gpt_agent import process_natural_language, generate_protocol, analyze_construct
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_process_natural_language():
    """Test natural language processing function"""
    logger.info("Testing natural language processing...")
    
    test_inputs = [
        "Insert GFP into pUC19 under T7 promoter",
        "Clone mCherry into pET28a with a lac promoter using Gibson Assembly",
        "Create a construct with RFP in pBluescript with trc promoter and add a 6xHis tag"
    ]
    
    for input_text in test_inputs:
        logger.info(f"Processing: '{input_text}'")
        result = process_natural_language(input_text)
        
        # Print the result
        logger.info(f"Result: {json.dumps(result, indent=2)}")
        
        # Basic validation
        assert "gene_name" in result, "Missing gene_name in result"
        assert "vector_name" in result, "Missing vector_name in result"
        assert "promoter_name" in result, "Missing promoter_name in result"
        
        logger.info(f"✅ Successfully processed: '{input_text}'")
    
    logger.info("Natural language processing test completed successfully!")

def test_generate_protocol():
    """Test protocol generation function"""
    logger.info("Testing protocol generation...")
    
    # Test with different construct configurations
    test_configs = [
        {
            "gene_name": "GFP",
            "vector_name": "pUC19",
            "promoter_name": "T7",
            "cloning_method": "Gibson Assembly"
        },
        {
            "gene_name": "mCherry",
            "vector_name": "pET28a",
            "promoter_name": "lac",
            "cloning_method": "Restriction Enzyme Cloning"
        }
    ]
    
    # Mock primer info
    primer_info = {
        "primers": [
            {
                "name": "Test_F",
                "sequence": "ATGCATGCATGCATGCATGC",
                "tm": 60.5,
                "gc_content": 50.0,
                "length": 20
            },
            {
                "name": "Test_R",
                "sequence": "GCATGCATGCATGCATGCAT",
                "tm": 60.5,
                "gc_content": 50.0,
                "length": 20
            }
        ]
    }
    
    for config in test_configs:
        logger.info(f"Generating protocol for: {config['gene_name']} in {config['vector_name']}")
        protocol = generate_protocol(config, primer_info)
        
        # Basic validation
        assert len(protocol) > 100, "Protocol is too short"
        assert "Materials" in protocol, "Protocol missing Materials section"
        assert "Procedure" in protocol, "Protocol missing Procedure section"
        
        # Write to file for inspection
        with open(f"test_protocol_{config['gene_name']}_{config['vector_name']}.md", "w") as f:
            f.write(protocol)
        
        logger.info(f"✅ Successfully generated protocol for: {config['gene_name']} in {config['vector_name']}")
    
    logger.info("Protocol generation test completed successfully!")

def test_analyze_construct():
    """Test construct analysis function"""
    logger.info("Testing construct analysis...")
    
    # Create a mock SeqRecord
    mock_record = SeqRecord(
        seq=Seq("ATGCATGCATGCATGCATGCATGCATGCATGC" * 100),
        id="GFP_pUC19",
        name="GFP_pUC19",
        description="GFP inserted into pUC19"
    )
    
    # Add annotations
    mock_record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    
    # Add simple features (would normally be more detailed)
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    mock_record.features = [
        SeqFeature(FeatureLocation(100, 900), type="CDS", qualifiers={"label": ["GFP"]}),
        SeqFeature(FeatureLocation(1000, 1500), type="promoter", qualifiers={"label": ["T7"]}),
        SeqFeature(FeatureLocation(2000, 2800), type="CDS", qualifiers={"label": ["Ampicillin resistance"]})
    ]
    
    # Test with different construct configurations
    test_configs = [
        {
            "gene_name": "GFP",
            "vector_name": "pUC19",
            "promoter_name": "T7",
            "cloning_method": "Gibson Assembly"
        },
        {
            "gene_name": "mCherry",
            "vector_name": "pET28a",
            "promoter_name": "lac",
            "cloning_method": "Restriction Enzyme Cloning"
        }
    ]
    
    for config in test_configs:
        logger.info(f"Analyzing construct: {config['gene_name']} in {config['vector_name']}")
        analysis = analyze_construct(mock_record, config)
        
        # Basic validation
        assert "potential_issues" in analysis, "Analysis missing potential_issues"
        assert "recommendations" in analysis, "Analysis missing recommendations"
        assert "expression_prediction" in analysis, "Analysis missing expression_prediction"
        
        # Print the result
        logger.info(f"Analysis result: {json.dumps(analysis, indent=2)}")
        
        logger.info(f"✅ Successfully analyzed construct: {config['gene_name']} in {config['vector_name']}")
    
    logger.info("Construct analysis test completed successfully!")

def run_all_tests():
    """Run all AI integration tests"""
    logger.info("Starting AI integration tests...")
    
    try:
        test_process_natural_language()
        test_generate_protocol()
        test_analyze_construct()
        
        logger.info("✅ All AI integration tests completed successfully!")
        return True
    except Exception as e:
        logger.error(f"❌ AI integration tests failed: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    # Check if API key is set
    api_key = os.environ.get("DEEPSEEK_API_KEY")
    if not api_key:
        logger.warning("No DeepSeek API key found in environment. Tests will use mock client.")
        logger.warning("Set the API key with: export DEEPSEEK_API_KEY=your_api_key")
    
    run_all_tests() 