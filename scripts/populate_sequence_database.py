import os
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import logging

# Add parent directory to path to import modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.sequence_database import NCBIRepository, AddgeneRepository, LocalRepository

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def fetch_ncbi_sequence(seq_id, output_path, format):
    # Implementation of fetch_ncbi_sequence function
    pass

def fetch_addgene_sequence(seq_id, output_path, format):
    # Implementation of fetch_addgene_sequence function
    pass

def fetch_igem_sequence(seq_id, output_path, format):
    # Implementation of fetch_igem_sequence function
    pass

def download_common_vectors():
    """Download common vectors from Addgene and NCBI"""
    # Create data directory if it doesn't exist
    data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
    os.makedirs(data_dir, exist_ok=True)
    
    # Initialize repositories
    ncbi_repo = NCBIRepository()
    addgene_repo = AddgeneRepository()
    local_repo = LocalRepository(data_dir)
    
    # Define vectors to download
    vectors = [
        {"id": "pUC19", "source": "ncbi", "accession": "L09137"},
        {"id": "pET28a", "source": "addgene", "accession": "69864"},
        {"id": "pET22b", "source": "addgene", "accession": "69744"},
        {"id": "pBluescript", "source": "ncbi", "accession": "X52327"},
    ]
    
    # Download each vector
    for vector in vectors:
        logger.info(f"Downloading {vector['id']} from {vector['source']}")
        
        try:
            # Fetch from the appropriate repository
            if vector['source'] == 'ncbi':
                record = ncbi_repo.fetch_sequence(vector['accession'])
            elif vector['source'] == 'addgene':
                record = addgene_repo.fetch_sequence(vector['accession'])
            else:
                logger.error(f"Unknown source: {vector['source']}")
                continue
            
            if not record:
                logger.error(f"Failed to fetch {vector['id']}")
                continue
            
            # Set the ID to the common name
            record.id = vector['id']
            record.name = vector['id']
            
            # Ensure annotations
            if not hasattr(record, 'annotations') or not record.annotations:
                record.annotations = {}
            record.annotations['molecule_type'] = 'DNA'
            record.annotations['topology'] = 'circular'
            
            # Save to local database
            filepath = os.path.join(data_dir, f"{vector['id']}.gb")
            SeqIO.write(record, filepath, "genbank")
            
            # Also save as FASTA
            fasta_path = os.path.join(data_dir, f"{vector['id']}.fasta")
            SeqIO.write(record, fasta_path, "fasta")
            
        except Exception as e:
            logger.error(f"Error downloading {vector['id']}: {str(e)}")

def main():
    # Promoters and regulatory elements
    sequences = [
        {"source": "ncbi", "id": "M13241.1", "name": "T7_promoter", "path": "data/promoters/T7_promoter.gb", "format": "genbank"},
        {"source": "ncbi", "id": "J01636.1", "name": "lac_promoter", "path": "data/promoters/lac_promoter.gb", "format": "genbank"},
        {"source": "igem", "id": "J56012", "name": "trc_promoter", "path": "data/promoters/trc_promoter.fasta", "format": "fasta"},
        {"source": "igem", "id": "I0500", "name": "araBAD_promoter", "path": "data/promoters/araBAD_promoter.fasta", "format": "fasta"},
        {"source": "igem", "id": "R0040", "name": "tet_promoter", "path": "data/promoters/tet_promoter.fasta", "format": "fasta"},
        {"source": "igem", "id": "B0015", "name": "terminator", "path": "data/terminators/B0015_terminator.fasta", "format": "fasta"},
        
        # Popular Addgene plasmids
        {"source": "addgene", "id": "42230", "name": "pX330", "path": "data/vectors/pX330.gb", "format": "genbank"},
        {"source": "addgene", "id": "48139", "name": "pSpCas9", "path": "data/vectors/pSpCas9.gb", "format": "genbank"},
        {"source": "addgene", "id": "17448", "name": "pLenti", "path": "data/vectors/pLenti.gb", "format": "genbank"},
        
        # Special parts
        {"source": "custom", "id": "6xHis", "name": "6xHis_tag", "path": "data/parts/6xHis_tag.fasta", "format": "fasta", 
         "sequence": "CATCATCATCATCATCAT", "description": "6xHis affinity tag"},
        {"source": "custom", "id": "FLAG", "name": "FLAG_tag", "path": "data/parts/FLAG_tag.fasta", "format": "fasta", 
         "sequence": "GACTACAAGGACGACGATGACAAG", "description": "FLAG affinity tag"},
        {"source": "custom", "id": "Strep", "name": "Strep_tag", "path": "data/parts/Strep_tag.fasta", "format": "fasta", 
         "sequence": "TGGTCTCATCCTCAATTTGAAAAG", "description": "Strep-tag II affinity tag"},
    ]
    
    # Track success/failure counts
    success_count = 0
    failure_count = 0
    
    # Process each sequence
    for seq_info in sequences:
        try:
            source = seq_info["source"]
            seq_id = seq_info["id"]
            output_path = seq_info["path"]
            format = seq_info.get("format", "genbank")
            
            # Skip if file already exists
            if os.path.exists(output_path):
                logger.info(f"Skipping {seq_info['name']}, file already exists: {output_path}")
                success_count += 1
                continue
            
            # Fetch from the appropriate source
            if source == "ncbi":
                success = fetch_ncbi_sequence(seq_id, output_path, format)
            elif source == "addgene":
                success = fetch_addgene_sequence(seq_id, output_path, format)
            elif source == "igem":
                success = fetch_igem_sequence(seq_id, output_path, format)
            elif source == "custom":
                # Create custom sequence
                record = SeqRecord(
                    Seq(seq_info["sequence"]),
                    id=seq_info["id"],
                    name=seq_info["name"],
                    description=seq_info["description"]
                )
                record.annotations = {'molecule_type': 'DNA'}
                
                # Ensure the output directory exists
                os.makedirs(os.path.dirname(output_path), exist_ok=True)
                
                # Save to file
                if format == "genbank":
                    SeqIO.write(record, output_path, "genbank")
                else:
                    SeqIO.write(record, output_path, "fasta")
                
                logger.info(f"Created custom sequence {seq_info['name']} at {output_path}")
                success = True
            else:
                logger.error(f"Unknown source: {source}")
                success = False
            
            # Track success/failure
            if success:
                success_count += 1
            else:
                failure_count += 1
                
            # Add a short delay between requests to avoid overwhelming APIs
            time.sleep(1)
            
        except Exception as e:
            logger.error(f"Error processing {seq_info.get('name', 'unknown')}: {str(e)}")
            failure_count += 1
    
    logger.info(f"Completed: {success_count} successes, {failure_count} failures")

if __name__ == "__main__":
    download_common_vectors()
    main() 