#!/usr/bin/env python3
"""
Test script for Ensilica to verify component functionality
"""
import os
import json
import logging
import argparse
from modules.gpt_agent import (
    process_natural_language, 
    generate_protocol, 
    test_api_connection, 
    DeepSeekClient
)
from modules.construct_builder import SequenceDatabase, ConstructBuilder
from modules.primer_tools import design_primers
from modules.visualize import visualize_construct
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def test_deepseek_connection():
    """Test the DeepSeek API connection"""
    logger.info("Testing DeepSeek API connection...")
    result = test_api_connection()
    if result:
        logger.info("✅ DeepSeek API connection successful")
    else:
        logger.error("❌ DeepSeek API connection failed")
    return result

def test_process_natural_language():
    """Test natural language processing"""
    logger.info("Testing natural language processing...")
    test_inputs = [
        "Insert GFP into pUC19 under a T7 promoter",
        "Clone mCherry into pET28a with a lac promoter"
    ]
    
    for input_text in test_inputs:
        logger.info(f"Processing: '{input_text}'")
        try:
            result = process_natural_language(input_text)
            logger.info(f"✅ Successfully processed: {json.dumps(result, indent=2)}")
        except Exception as e:
            logger.error(f"❌ Failed to process: {str(e)}")
            return False
    
    return True

def test_sequence_database():
    """Test sequence database functionality"""
    logger.info("Testing sequence database...")
    db = SequenceDatabase()
    
    # Test loading sequences
    sequences_to_test = [
        ("pUC19", "vector"),
        ("GFP", "gene"),
        ("mCherry", "gene")
    ]
    
    for seq_name, seq_type in sequences_to_test:
        logger.info(f"Loading {seq_type}: {seq_name}")
        try:
            seq_record = db.get_sequence(seq_name, seq_type)
            logger.info(f"✅ Successfully loaded {seq_name} ({len(seq_record.seq)} bp)")
        except Exception as e:
            logger.error(f"❌ Failed to load {seq_name}: {str(e)}")
            return False
    
    return True

def test_construct_builder():
    """Test construct builder functionality"""
    logger.info("Testing construct builder...")
    db = SequenceDatabase()
    builder = ConstructBuilder(db)
    
    # Test construct creation
    test_constructs = [
        {
            "gene_name": "GFP",
            "vector_name": "pUC19",
            "promoter_name": "T7",
            "cloning_method": "Gibson Assembly"
        },
        {
            "gene_name": "mCherry",
            "vector_name": "pUC19",
            "promoter_name": "lac",
            "cloning_method": "Restriction Enzyme Cloning"
        }
    ]
    
    for construct_info in test_constructs:
        construct_desc = f"{construct_info['gene_name']} in {construct_info['vector_name']}"
        logger.info(f"Building construct: {construct_desc}")
        try:
            result = builder.create_construct(construct_info)
            logger.info(f"✅ Successfully built {construct_desc} ({len(result['construct_record'].seq)} bp)")
        except Exception as e:
            logger.error(f"❌ Failed to build {construct_desc}: {str(e)}")
            return False
    
    return True

def test_primer_design():
    """Test primer design functionality"""
    logger.info("Testing primer design...")
    
    # Mock insert region
    insert_region = {
        "start": 100,
        "end": 400,
        "sequence": "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCC"
    }
    
    # Mock construct info
    construct_info = {
        "gene_name": "GFP",
        "vector_name": "pUC19",
        "promoter_name": "T7",
        "cloning_method": "Gibson Assembly"
    }
    
    try:
        logger.info(f"Designing primers for {construct_info['gene_name']}")
        result = design_primers(insert_region, construct_info)
        logger.info(f"✅ Successfully designed {len(result['primers'])} primers")
        
        # Log the primers
        for primer in result['primers']:
            logger.info(f"  - {primer['name']}: {primer['sequence']} (Tm: {primer['tm']}°C)")
            
    except Exception as e:
        logger.error(f"❌ Failed to design primers: {str(e)}")
        return False
    
    return True

def test_visualization():
    """Test visualization functionality"""
    logger.info("Testing visualization...")
    
    # Create a simple mock construct
    mock_seq = Seq("ATGC" * 500)  # 2000 bp
    mock_record = SeqRecord(mock_seq, id="test_construct", name="test_construct")
    
    # Mock insert region
    mock_insert = {
        "start": 500,
        "end": 800,
        "sequence": "ATGC" * 75
    }
    
    try:
        logger.info("Generating visualization")
        result = visualize_construct(mock_record, mock_insert)
        logger.info(f"✅ Successfully generated visualization ({len(result)} characters)")
    except Exception as e:
        logger.error(f"❌ Failed to generate visualization: {str(e)}")
        return False
    
    return True

def test_protocol_generation():
    """Test protocol generation"""
    logger.info("Testing protocol generation...")
    
    # Mock construct info
    construct_info = {
        "gene_name": "GFP",
        "vector_name": "pUC19",
        "promoter_name": "T7",
        "cloning_method": "Gibson Assembly"
    }
    
    # Mock primer info
    primer_info = {
        "primers": [
            {
                "name": "GFP_F",
                "sequence": "ATGGTGAGCAAGGGCGAGGAG",
                "tm": 64.5,
                "gc_content": 57.1,
                "length": 21
            },
            {
                "name": "GFP_R",
                "sequence": "TTACTTGTACAGCTCGTCCAT",
                "tm": 58.3,
                "gc_content": 42.9,
                "length": 21
            }
        ]
    }
    
    try:
        logger.info(f"Generating protocol for {construct_info['gene_name']} in {construct_info['vector_name']}")
        result = generate_protocol(construct_info, primer_info)
        logger.info(f"✅ Successfully generated protocol ({len(result)} characters)")
    except Exception as e:
        logger.error(f"❌ Failed to generate protocol: {str(e)}")
        return False
    
    return True

def run_all_tests():
    """Run all tests and report results"""
    tests = [
        ("DeepSeek API Connection", test_deepseek_connection),
        ("Natural Language Processing", test_process_natural_language),
        ("Sequence Database", test_sequence_database),
        ("Construct Builder", test_construct_builder),
        ("Primer Design", test_primer_design),
        ("Visualization", test_visualization),
        ("Protocol Generation", test_protocol_generation)
    ]
    
    results = {}
    all_passed = True
    
    for test_name, test_func in tests:
        logger.info(f"\n{'='*50}\nRunning test: {test_name}\n{'='*50}")
        try:
            result = test_func()
            results[test_name] = result
            if not result:
                all_passed = False
        except Exception as e:
            logger.error(f"Test '{test_name}' failed with exception: {str(e)}")
            results[test_name] = False
            all_passed = False
    
    # Print summary
    logger.info("\n\n")
    logger.info("="*50)
    logger.info("TEST RESULTS SUMMARY")
    logger.info("="*50)
    
    for test_name, result in results.items():
        status = "✅ PASSED" if result else "❌ FAILED"
        logger.info(f"{test_name}: {status}")
    
    overall = "✅ ALL TESTS PASSED" if all_passed else "❌ SOME TESTS FAILED"
    logger.info(f"\nOverall: {overall}")
    
    return all_passed

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test Ensilica functionality")
    parser.add_argument("--test", choices=["all", "api", "nlp", "db", "builder", "primers", "viz", "protocol"], 
                      default="all", help="Specific test to run (default: all)")
    args = parser.parse_args()
    
    # Check API key
    if not os.environ.get("DEEPSEEK_API_KEY"):
        logger.warning("\n⚠️  DEEPSEEK_API_KEY environment variable is not set!")
        logger.warning("API-dependent tests will fail. Set the key with:")
        logger.warning("export DEEPSEEK_API_KEY=your_api_key_here\n")
    
    # Run requested test
    if args.test == "all":
        run_all_tests()
    elif args.test == "api":
        test_deepseek_connection()
    elif args.test == "nlp":
        test_process_natural_language()
    elif args.test == "db":
        test_sequence_database()
    elif args.test == "builder":
        test_construct_builder()
    elif args.test == "primers":
        test_primer_design()
    elif args.test == "viz":
        test_visualization()
    elif args.test == "protocol":
        test_protocol_generation() 