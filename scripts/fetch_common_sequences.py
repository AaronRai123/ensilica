#!/usr/bin/env python3
"""
Script to pre-populate the data directory with common sequences
"""
import os
import sys
import logging

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.sequence_database import SequenceDatabase

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    """Fetch common sequences and save to local repository"""
    # Ensure data directories exist
    os.makedirs("data/vectors", exist_ok=True)
    os.makedirs("data/genes", exist_ok=True)
    os.makedirs("data/parts", exist_ok=True)
    
    # Initialize sequence database
    db = SequenceDatabase()
    
    # List of sequences to fetch
    sequences = [
        # Vectors
        {"id": "pUC19", "repo": "ncbi", "alt_id": "M77789.2"},
        {"id": "pET28a", "repo": "ncbi", "alt_id": "JX462669.1"},
        {"id": "pBluescript", "repo": "ncbi", "alt_id": "X52329.1"},
        {"id": "pBAD", "repo": "ncbi", "alt_id": "X17768.1"},
        {"id": "pACYC184", "repo": "ncbi", "alt_id": "X05858.1"},
        {"id": "pBR322", "repo": "ncbi", "alt_id": "X03469.1"},
        {"id": "pRSF1b", "repo": "ncbi", "alt_id": "X03470.1"},
        {"id": "pTrc99a", "repo": "ncbi", "alt_id": "X03471.1"},
        {"id": "pBBR1MCS", "repo": "ncbi", "alt_id": "X03472.1"},
        {"id": "pSB1C3", "repo": "ncbi", "alt_id": "X03473.1"},
        # Genes and parts
        {"id": "GFP", "repo": "local", "alt_id": "GFP"},
        {"id": "mCherry", "repo": "local", "alt_id": "mCherry"},
        {"id": "RFP", "repo": "local", "alt_id": "RFP"},
        {"id": "YFP", "repo": "local", "alt_id": "YFP"},
        {"id": "BFP", "repo": "local", "alt_id": "BFP"},
        {"id": "kanR", "repo": "local", "alt_id": "kanR"},
        {"id": "ampR", "repo": "local", "alt_id": "ampR"},
        {"id": "camR", "repo": "local", "alt_id": "camR"},
        {"id": "tetR", "repo": "local", "alt_id": "tetR"},
        {"id": "T7", "repo": "local", "alt_id": "T7"},
        {"id": "lac", "repo": "local", "alt_id": "lac"},
        {"id": "trc", "repo": "local", "alt_id": "trc"},
        {"id": "araBAD", "repo": "local", "alt_id": "araBAD"},
        {"id": "tet", "repo": "local", "alt_id": "tet"},
        {"id": "6xHis", "repo": "local", "alt_id": "6xHis"},
        {"id": "FLAG", "repo": "local", "alt_id": "FLAG"},
        {"id": "MBP", "repo": "local", "alt_id": "MBP"},
        {"id": "GST", "repo": "local", "alt_id": "GST"},
        {"id": "Strep", "repo": "local", "alt_id": "Strep"},
        {"id": "lacZ", "repo": "local", "alt_id": "lacZ"},
        {"id": "GUS", "repo": "local", "alt_id": "GUS"},
        {"id": "luciferase", "repo": "local", "alt_id": "luciferase"}
    ]
    
    for seq in sequences:
        logger.info(f"Fetching sequence: {seq['id']}")
        record = db.get_sequence(seq['id'], seq['repo'])
        if record:
            logger.info(f"Successfully fetched {seq['id']} from {seq['repo']}")
            if 'alt_id' in seq:
                logger.info(f"Alternative ID: {seq['alt_id']}")
            db.repositories["local"].save_sequence(record)
        else:
            logger.warning(f"Failed to fetch {seq['id']} from {seq['repo']}")

if __name__ == "__main__":
    main() 