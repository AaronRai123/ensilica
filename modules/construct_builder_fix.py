"""
Construct Builder module for creating DNA constructs from structured information
with fixes for Biopython GenBank export
"""
from typing import Dict, Any
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import io

# Define paths to sequence data files
DATA_DIR = "data"

def load_sequence(sequence_name: str, sequence_type: str) -> SeqRecord: 