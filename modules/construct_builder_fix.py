"""
Construct Builder module for creating DNA constructs from structured information
with fixes for Biopython GenBank export
"""
from typing import Dict, Any, Tuple, List, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import io
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define paths to sequence data files
DATA_DIR = "data"

# Export these functions for use in other modules
__all__ = ['load_sequence', 'ensure_valid_annotations', 'create_construct']

def load_sequence(sequence_name: str, sequence_type: str) -> SeqRecord:
    """
    Load a sequence from a local file based on name and type
    
    Args:
        sequence_name: Name of the sequence
        sequence_type: Type of sequence (gene, vector, etc.)
        
    Returns:
        SeqRecord object
    """
    # Determine file path based on sequence type
    if sequence_type == "gene":
        subdir = "genes"
    elif sequence_type == "vector":
        subdir = "vectors"
    else:
        subdir = ""
    
    # Check different possible file locations
    possible_paths = [
        os.path.join(DATA_DIR, subdir, f"{sequence_name}.gb"),
        os.path.join(DATA_DIR, subdir, f"{sequence_name}.fasta"),
        os.path.join(DATA_DIR, f"{sequence_name}.gb"),
        os.path.join(DATA_DIR, f"{sequence_name}.fasta"),
    ]
    
    # Try each path
    for path in possible_paths:
        if os.path.exists(path):
            # Determine format from extension
            file_format = "genbank" if path.endswith(".gb") else "fasta"
            record = SeqIO.read(path, file_format)
            
            # Ensure proper annotations for GenBank export
            ensure_valid_annotations(record)
            
            return record
    
    # If sequence not found
    raise FileNotFoundError(f"Sequence {sequence_name} not found in {DATA_DIR}")

def ensure_valid_annotations(record: SeqRecord) -> None:
    """
    Ensure the SeqRecord has valid annotations for GenBank export
    
    Args:
        record: SeqRecord to validate/fix
    """
    # Initialize annotations dict if needed
    if not hasattr(record, 'annotations') or record.annotations is None:
        record.annotations = {}
    
    # Ensure molecule_type is set (required for GenBank export)
    if 'molecule_type' not in record.annotations:
        record.annotations['molecule_type'] = 'DNA'
    
    # Set topology for plasmids if not specified
    if 'topology' not in record.annotations:
        record.annotations['topology'] = 'circular'
    
    # Add date if missing
    if 'date' not in record.annotations:
        from datetime import datetime
        record.annotations['date'] = datetime.now().strftime("%d-%b-%Y").upper()

def create_construct(construct_info: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create a DNA construct from structured information
    
    Args:
        construct_info: Dictionary with construct information including:
            - gene_name: Name of the gene
            - vector_name: Name of the vector
            - cloning_method: Method for cloning 
            - promoter_name: Optional promoter name
            - insertion_site: Optional insertion site
            
    Returns:
        Dictionary with construct information
    """
    # Extract required info
    gene_name = construct_info.get("gene_name")
    vector_name = construct_info.get("vector_name")
    cloning_method = construct_info.get("cloning_method", "Gibson Assembly")
    
    if not gene_name or not vector_name:
        raise ValueError("Missing required gene_name or vector_name in construct_info")
    
    # Load sequences
    try:
        gene_record = load_sequence(gene_name, "gene")
        vector_record = load_sequence(vector_name, "vector")
        
        # Create a new construct record
        construct_id = f"{gene_name}_{vector_name}_{cloning_method.lower().replace(' ', '_')}"
        construct_record = SeqRecord(
            Seq(str(vector_record.seq) + str(gene_record.seq)),
            id=construct_id,
            name=construct_id,
            description=f"Construct with {gene_name} in {vector_name} using {cloning_method}"
        )
        
        # Ensure proper annotations for GenBank export
        ensure_valid_annotations(construct_record)
        
        # Add features from both records
        for feature in vector_record.features + gene_record.features:
            if feature.type != "source":
                construct_record.features.append(feature)
        
        # Add a source feature for the entire construct
        source_feature = SeqFeature(
            FeatureLocation(0, len(construct_record.seq)),
            type="source",
            qualifiers={"organism": "synthetic construct"}
        )
        construct_record.features.insert(0, source_feature)
        
        # Generate GenBank and FASTA content
        gb_handle = io.StringIO()
        SeqIO.write(construct_record, gb_handle, "genbank")
        genbank_content = gb_handle.getvalue()
        
        fasta_handle = io.StringIO()
        SeqIO.write(construct_record, fasta_handle, "fasta")
        fasta_content = fasta_handle.getvalue()
        
        # Extract the insert region for primer design
        insert_region = {
            "sequence": str(gene_record.seq),
            "start": len(vector_record.seq),
            "end": len(vector_record.seq) + len(gene_record.seq)
        }
        
        return {
            "construct_record": construct_record,
            "genbank_content": genbank_content,
            "fasta_content": fasta_content,
            "gene_info": {"name": gene_name},
            "vector_info": {"name": vector_name},
            "cloning_method": cloning_method,
            "insert_region": insert_region
        }
    
    except Exception as e:
        logger.error(f"Error creating construct: {str(e)}")
        raise RuntimeError(f"Failed to create construct: {str(e)}") from e 