"""
Enhanced Construct Builder module with support for real sequences and advanced cloning methods
"""
from typing import Dict, Any, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import io
import logging
from modules.sequence_database import get_sequence, search_sequences, get_available_sequences
from modules.cloning_methods import (
    gibson_assembly, 
    restriction_cloning, 
    golden_gate_assembly, 
    gateway_cloning,
    slic_assembly,
    find_restriction_sites
)

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ConstructBuilder:
    """Builder for DNA constructs using various cloning methods"""
    
    def __init__(self, data_dir: str = "data"):
        """
        Initialize the construct builder
        
        Args:
            data_dir: Directory containing sequence data files
        """
        self.data_dir = data_dir
    
    def create_construct(self, construct_info: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create a DNA construct based on the structured information
        
        Args:
            construct_info: Dictionary containing structured information about the construct
            
        Returns:
            Dictionary containing the constructed sequence and related information
        
        Raises:
            ValueError: If required sequences can't be found or are invalid
        """
        # Extract key information
        gene_id = construct_info.get("gene_name", "unknown_gene")
        vector_id = construct_info.get("vector_name", "pUC19")
        promoter_id = construct_info.get("promoter_name", "lac")
        cloning_method = construct_info.get("cloning_method", "Gibson Assembly")
        insertion_site = construct_info.get("insertion_site", "Multiple Cloning Site")
        additional_features = construct_info.get("additional_features", [])
        
        # Ensure additional_features is a list
        if additional_features is None:
            additional_features = []
        
        try:
            # Load the sequences
            logger.info(f"Loading gene: {gene_id}")
            gene = self._load_sequence(gene_id, "gene")
            
            logger.info(f"Loading vector: {vector_id}")
            vector = self._load_sequence(vector_id, "vector")
            
            # Create the construct using the appropriate cloning method
            logger.info(f"Creating construct using {cloning_method}")
            
            if cloning_method == "Gibson Assembly":
                construct_record, metadata = gibson_assembly(vector, gene)
            elif cloning_method == "Restriction Enzyme Cloning":
                # For restriction cloning, identify common enzymes that don't cut the insert
                gene_sites = find_restriction_sites(gene)
                vector_sites = find_restriction_sites(vector)
                
                # Find enzymes that cut the vector but not the gene
                suitable_enzymes = []
                for enzyme, sites in vector_sites.items():
                    if enzyme not in gene_sites or not gene_sites[enzyme]:
                        if sites:  # Ensure the enzyme cuts the vector
                            suitable_enzymes.append(enzyme)
                
                # Use the first two suitable enzymes, or fall back to EcoRI/HindIII
                if len(suitable_enzymes) >= 2:
                    enzymes = suitable_enzymes[:2]
                else:
                    enzymes = ['EcoRI', 'HindIII']
                    logger.warning(f"No ideal restriction enzymes found, using default: {', '.join(enzymes)}")
                
                construct_record, metadata = restriction_cloning(vector, gene, enzymes)
            elif cloning_method == "Golden Gate":
                # For Golden Gate, we treat the gene as a single insert
                construct_record, metadata = golden_gate_assembly(vector, [gene])
            elif cloning_method == "Gateway Cloning":
                # For Gateway, we would need proper entry and destination vectors
                construct_record, metadata = gateway_cloning(gene, vector)
            elif cloning_method == "SLIC":
                construct_record, metadata = slic_assembly(vector, gene)
            else:
                # Default to Gibson Assembly
                logger.warning(f"Unknown cloning method: {cloning_method}, falling back to Gibson Assembly")
                construct_record, metadata = gibson_assembly(vector, gene)
            
            # Add promoter if specified
            if promoter_id:
                self._add_promoter(construct_record, promoter_id, metadata.get("insertion_site", 0))
            
            # Add any additional features
            for feature_info in additional_features:
                self._add_feature(construct_record, feature_info)
            
            # Convert to different formats
            gb_output = io.StringIO()
            SeqIO.write(construct_record, gb_output, "genbank")
            
            fasta_output = io.StringIO()
            SeqIO.write(construct_record, fasta_output, "fasta")
            
            # Prepare result
            result = {
                "construct_record": construct_record,
                "genbank_content": gb_output.getvalue(),
                "fasta_content": fasta_output.getvalue(),
                "insert_region": {
                    "start": metadata.get("insertion_site", 0),
                    "end": metadata.get("insertion_site", 0) + len(gene.seq),
                    "sequence": str(gene.seq)
                },
                "cloning_metadata": metadata
            }
            
            logger.info(f"Successfully created {construct_record.id} ({len(construct_record.seq)} bp)")
            return result
            
        except Exception as e:
            error_msg = f"Failed to create construct: {str(e)}"
            logger.error(error_msg)
            raise ValueError(error_msg) from e
    
    def _load_sequence(self, sequence_id: str, sequence_type: str) -> SeqRecord:
        """
        Load a sequence by ID from the sequence database
        
        Args:
            sequence_id: ID of the sequence
            sequence_type: Type of sequence (gene, vector, etc.)
            
        Returns:
            SeqRecord object
        
        Raises:
            ValueError: If sequence is not found or invalid
        """
        # Try to get sequence from local data first
        record = get_sequence(sequence_id)
        
        if record:
            logger.info(f"Using local sequence for {sequence_id}")
            # Validate the sequence
            if sequence_type == 'gene':
                self._validate_sequence(record, sequence_id)
            
            # Ensure required annotations
            if not hasattr(record, 'annotations') or not record.annotations:
                record.annotations = {}
            if 'molecule_type' not in record.annotations:
                record.annotations['molecule_type'] = 'DNA'
            if sequence_type == 'vector' and 'topology' not in record.annotations:
                record.annotations['topology'] = 'circular'
            
            return record
        
        # If not found locally, try external repositories
        if sequence_type == 'vector':
            # Try to fetch common vectors from Addgene or NCBI
            logger.info(f"Trying to fetch vector {sequence_id} from external repositories")
            
            # Common vector prefixes
            if sequence_id.startswith("pET"):
                # Try Addgene for pET vectors
                from modules.sequence_database import AddgeneRepository
                repo = AddgeneRepository()
                addgene_id = None
                
                # Map common pET vectors to Addgene IDs
                if sequence_id == "pET28a":
                    addgene_id = "69864"  # pET28a(+)
                elif sequence_id == "pET22b":
                    addgene_id = "69744"  # pET22b(+)
                
                if addgene_id:
                    record = repo.fetch_sequence(addgene_id)
                    if record:
                        logger.info(f"Successfully fetched {sequence_id} from Addgene (ID: {addgene_id})")
                        # Update ID to match requested ID
                        record.id = sequence_id
                        record.name = sequence_id
                        
                        if sequence_type == 'gene':
                            self._validate_sequence(record, sequence_id)
                        
                        return record
            
            elif sequence_id == "pUC19":
                # Try NCBI for pUC19
                from modules.sequence_database import NCBIRepository
                repo = NCBIRepository()
                record = repo.fetch_sequence("L09137")  # NCBI accession for pUC19
                if record:
                    logger.info(f"Successfully fetched {sequence_id} from NCBI")
                    # Update ID to match requested ID
                    record.id = sequence_id
                    record.name = sequence_id
                    return record
        
        elif sequence_type == 'gene':
            # Try to fetch genes from NCBI
            logger.info(f"Trying to fetch gene {sequence_id} from external repositories")
            
            from modules.sequence_database import NCBIRepository
            repo = NCBIRepository()
            
            # Try direct accession lookup for common genes
            accession = None
            if sequence_id == "GFP":
                accession = "U55762"  # GFP from Aequorea victoria
            elif sequence_id == "mCherry":
                accession = "AY678264"  # mCherry red fluorescent protein
            elif sequence_id == "RFP":
                accession = "AF506027"  # DsRed-Express red fluorescent protein
            elif sequence_id == "TFEB":
                accession = "NM_001167827"  # Human TFEB
            
            if accession:
                record = repo.fetch_sequence(accession)
                if record:
                    logger.info(f"Successfully fetched {sequence_id} from NCBI (accession: {accession})")
                    # Update ID to match requested ID
                    record.id = sequence_id
                    record.name = sequence_id
                    
                    # Extract the CDS feature if it exists
                    cds_features = [f for f in record.features if f.type == "CDS"]
                    if cds_features:
                        # Extract the CDS sequence
                        cds_feature = cds_features[0]
                        start = cds_feature.location.start
                        end = cds_feature.location.end
                        cds_seq = record.seq[start:end]
                        
                        # Create a new record with just the CDS
                        from Bio.SeqRecord import SeqRecord
                        from Bio.Seq import Seq
                        cds_record = SeqRecord(
                            cds_seq,
                            id=record.id,
                            name=record.name,
                            description=f"{sequence_id} coding sequence"
                        )
                        cds_record.annotations = {'molecule_type': 'DNA'}
                        
                        # Add CDS feature
                        from Bio.SeqFeature import SeqFeature, FeatureLocation
                        cds_feature = SeqFeature(
                            FeatureLocation(0, len(cds_seq), 1),
                            type="CDS",
                            qualifiers={"label": [sequence_id], "gene": [sequence_id]}
                        )
                        cds_record.features = [cds_feature]
                        
                        # Validate the sequence
                        self._validate_sequence(cds_record, sequence_id)
                        
                        return cds_record
                    
                    # Validate the sequence
                    self._validate_sequence(record, sequence_id)
                    return record
            
            # If not found by direct accession, try search
            logger.info(f"No direct accession for {sequence_id}, attempting search")
            search_results = search_sequences(f"{sequence_id}[Gene Name]", ["NCBI"], max_results=1)
            
            if search_results:
                search_hit = search_results[0]
                record_id = search_hit.get("id")
                
                if record_id:
                    record = repo.fetch_sequence(record_id)
                    if record:
                        logger.info(f"Successfully fetched {sequence_id} from NCBI search (ID: {record_id})")
                        # Update ID to match requested ID
                        record.id = sequence_id
                        record.name = sequence_id
                        
                        # Extract CDS if possible
                        cds_features = [f for f in record.features if f.type == "CDS"]
                        if cds_features:
                            cds_feature = cds_features[0]
                            start = cds_feature.location.start
                            end = cds_feature.location.end
                            cds_seq = record.seq[start:end]
                            
                            cds_record = SeqRecord(
                                cds_seq,
                                id=record.id,
                                name=record.name,
                                description=f"{sequence_id} coding sequence"
                            )
                            cds_record.annotations = {'molecule_type': 'DNA'}
                            
                            # Add CDS feature
                            cds_feature = SeqFeature(
                                FeatureLocation(0, len(cds_seq), 1),
                                type="CDS",
                                qualifiers={"label": [sequence_id], "gene": [sequence_id]}
                            )
                            cds_record.features = [cds_feature]
                            
                            # Validate the sequence
                            self._validate_sequence(cds_record, sequence_id)
                            
                            return cds_record
                        
                        # Validate the sequence
                        self._validate_sequence(record, sequence_id)
                        return record
        
        # If we get here, we couldn't find the sequence
        error_msg = f"Could not find sequence for {sequence_id} in any repository."
        logger.error(error_msg)
        raise ValueError(error_msg)
    
    def _validate_sequence(self, record: SeqRecord, sequence_id: str) -> None:
        """
        Validate a gene sequence to ensure it meets quality requirements
        
        Args:
            record: SeqRecord object to validate
            sequence_id: ID of the sequence for error messages
            
        Raises:
            ValueError: If sequence is invalid
        """
        if not record or not hasattr(record, 'seq') or not record.seq:
            raise ValueError(f"Invalid sequence for {sequence_id}: empty sequence")
        
        # Check sequence length
        if len(record.seq) < 60:
            raise ValueError(f"Invalid sequence for {sequence_id}: too short (minimum 60 bp required)")
        
        # Count ambiguous bases (N)
        seq_str = str(record.seq).upper()
        n_count = seq_str.count('N')
        n_percentage = (n_count / len(seq_str)) * 100
        
        if n_percentage > 10:
            raise ValueError(f"Invalid sequence for {sequence_id}: contains {n_percentage:.1f}% ambiguous bases (maximum 10% allowed)")
        
        # Smart sequence validation
        validation_results = self.smart_sequence_validation(record)
        
        # Log the validation results
        for key, result in validation_results.items():
            if not result["valid"]:
                logger.warning(f"Sequence {sequence_id} validation issue: {result['message']}")
        
        logger.info(f"Validated sequence for {sequence_id}: {len(record.seq)} bp, {n_percentage:.1f}% ambiguous bases")
    
    def smart_sequence_validation(self, record: SeqRecord) -> Dict[str, Any]:
        """
        Perform smart sequence validation checks on a SeqRecord
        
        Args:
            record: SeqRecord object to validate
            
        Returns:
            Dictionary with validation results for different aspects
        """
        results = {
            "frame_shift": {"valid": True, "message": ""},
            "promoter_orientation": {"valid": True, "message": ""},
            "rbs_presence": {"valid": True, "message": ""},
            "restriction_sites": {"valid": True, "message": "", "conflicting_sites": []}
        }
        
        seq_str = str(record.seq).upper()
        
        # 1. Check for frame shift (start/stop codon)
        results["frame_shift"] = self._check_frame_shift(record)
        
        # 2. Check promoter orientation
        results["promoter_orientation"] = self._check_promoter_orientation(record)
        
        # 3. Check for RBS presence before CDS
        results["rbs_presence"] = self._check_rbs_presence(record)
        
        # 4. Check for internal restriction sites that conflict with common cloning enzymes
        results["restriction_sites"] = self._check_restriction_sites(record)
        
        return results
        
    def _check_frame_shift(self, record: SeqRecord) -> Dict[str, Any]:
        """Check for proper reading frame with start/stop codons"""
        result = {"valid": True, "message": ""}
        seq_str = str(record.seq).upper()
        
        # Look for CDS features first
        cds_features = [f for f in record.features if f.type == "CDS"]
        
        if cds_features:
            for i, feature in enumerate(cds_features):
                # Extract the feature sequence
                feature_seq = feature.extract(record.seq)
                feature_str = str(feature_seq).upper()
                
                # Check for start codon
                if not feature_str.startswith('ATG'):
                    result["valid"] = False
                    result["message"] += f"CDS #{i+1} does not start with ATG. "
                
                # Check for in-frame stop codon
                has_stop = False
                for j in range(0, len(feature_str) - 2, 3):
                    codon = feature_str[j:j+3]
                    if codon in ['TAA', 'TAG', 'TGA']:
                        if j == len(feature_str) - 3:  # Last codon
                            has_stop = True
                        else:
                            result["valid"] = False
                            result["message"] += f"CDS #{i+1} has premature stop codon at position {j}. "
                
                if not has_stop:
                    result["valid"] = False
                    result["message"] += f"CDS #{i+1} has no stop codon. "
        else:
            # If no CDS features, do a basic check on the whole sequence
            has_start = 'ATG' in seq_str
            has_stop = any(stop in seq_str for stop in ['TAA', 'TAG', 'TGA'])
            
            if not has_start:
                result["valid"] = False
                result["message"] += "No start codon (ATG) found. "
            
            if not has_stop:
                result["valid"] = False
                result["message"] += "No stop codon (TAA/TAG/TGA) found. "
        
        return result
    
    def _check_promoter_orientation(self, record: SeqRecord) -> Dict[str, Any]:
        """Check for proper promoter orientation relative to CDS"""
        result = {"valid": True, "message": ""}
        
        promoter_features = [f for f in record.features if f.type == "promoter"]
        cds_features = [f for f in record.features if f.type == "CDS"]
        
        if promoter_features and cds_features:
            for promoter in promoter_features:
                promoter_end = promoter.location.end
                
                # Find CDS features that should be downstream of this promoter
                for i, cds in enumerate(cds_features):
                    cds_start = cds.location.start
                    
                    # If the CDS is on reverse strand but promoter isn't, or vice versa
                    if promoter.strand != cds.strand:
                        result["valid"] = False
                        result["message"] += f"Promoter and CDS #{i+1} have incompatible orientations (different strands). "
                    
                    # Check if the promoter is upstream of the CDS (within reasonable distance)
                    if promoter.strand == 1:  # Forward strand
                        if not (promoter_end <= cds_start and cds_start - promoter_end < 500):
                            result["valid"] = False
                            result["message"] += f"Promoter not properly positioned upstream of CDS #{i+1}. "
                    else:  # Reverse strand
                        if not (cds_start <= promoter_end and promoter_end - cds_start < 500):
                            result["valid"] = False
                            result["message"] += f"Promoter not properly positioned upstream of CDS #{i+1} on reverse strand. "
        
        return result
    
    def _check_rbs_presence(self, record: SeqRecord) -> Dict[str, Any]:
        """Check for RBS presence before CDS"""
        result = {"valid": True, "message": ""}
        
        cds_features = [f for f in record.features if f.type == "CDS"]
        rbs_features = [f for f in record.features if f.type in ["RBS", "ribosome_binding_site"]]
        
        # Common RBS consensus sequences to check for
        rbs_consensus = ["AGGAGG", "GGAGG", "AGGAG", "GGAG", "GAGG", "AGGA", "AAGG"]
        
        if cds_features:
            for i, cds in enumerate(cds_features):
                has_rbs = False
                
                # First check if there's an RBS feature associated with this CDS
                for rbs in rbs_features:
                    if cds.strand == 1:  # Forward strand
                        if 0 < cds.location.start - rbs.location.end < 20:
                            has_rbs = True
                            break
                    else:  # Reverse strand
                        if 0 < rbs.location.start - cds.location.end < 20:
                            has_rbs = True
                            break
                
                # If no RBS feature, look for RBS sequences in the upstream region
                if not has_rbs:
                    if cds.strand == 1:  # Forward strand
                        upstream_start = max(0, cds.location.start - 20)
                        upstream_end = cds.location.start
                        upstream_seq = str(record.seq[upstream_start:upstream_end]).upper()
                        
                        for consensus in rbs_consensus:
                            if consensus in upstream_seq:
                                has_rbs = True
                                break
                    else:  # Reverse strand
                        upstream_start = cds.location.end
                        upstream_end = min(len(record.seq), cds.location.end + 20)
                        upstream_seq = str(record.seq[upstream_start:upstream_end]).upper()
                        # Need reverse complement for the consensus on reverse strand
                        from Bio.Seq import Seq
                        for consensus in rbs_consensus:
                            rev_consensus = str(Seq(consensus).reverse_complement())
                            if rev_consensus in upstream_seq:
                                has_rbs = True
                                break
                
                if not has_rbs:
                    result["valid"] = False
                    result["message"] += f"No RBS found upstream of CDS #{i+1}. "
        
        return result
    
    def _check_restriction_sites(self, record: SeqRecord) -> Dict[str, Any]:
        """Check for restriction sites that could interfere with cloning"""
        result = {"valid": True, "message": "", "conflicting_sites": []}
        
        # Common restriction enzymes used for cloning
        common_enzymes = {
            "EcoRI": "GAATTC",
            "BamHI": "GGATCC",
            "HindIII": "AAGCTT",
            "XhoI": "CTCGAG",
            "XbaI": "TCTAGA",
            "PstI": "CTGCAG",
            "SalI": "GTCGAC",
            "SmaI": "CCCGGG",
            "KpnI": "GGTACC",
            "NdeI": "CATATG",
            "NotI": "GCGGCCGC",
            "BglII": "AGATCT"
        }
        
        seq_str = str(record.seq).upper()
        
        for enzyme, site in common_enzymes.items():
            # Count occurrences of the site in the sequence
            count = seq_str.count(site)
            
            if count > 0:
                result["conflicting_sites"].append({
                    "enzyme": enzyme,
                    "site": site,
                    "count": count,
                    "positions": [seq_str.find(site, i) for i in range(0, len(seq_str), len(site)) if seq_str.find(site, i) != -1]
                })
        
        # If we found restriction sites that might interfere with cloning
        if result["conflicting_sites"]:
            # Still valid but with a warning
            result["message"] = f"Found {len(result['conflicting_sites'])} restriction sites that might interfere with cloning: "
            result["message"] += ", ".join([f"{site['enzyme']} ({site['count']})" for site in result["conflicting_sites"]])
        
        return result
    
    def _create_mock_sequence(self, sequence_id: str, sequence_type: str) -> SeqRecord:
        """
        Create a mock sequence when a real one is not available
        
        Args:
            sequence_id: ID for the mock sequence
            sequence_type: Type of sequence to create
            
        Returns:
            SeqRecord object with mock sequence
        """
        logger.info(f"Creating mock {sequence_type} sequence for {sequence_id}")
        
        if sequence_type == "vector":
            # Create a mock vector sequence (about 3kb)
            seq = Seq("ATGCATGCATGCATGCATGCATGCATGCATGCATGC" * 80)  # ~3kb
            record = SeqRecord(seq, id=sequence_id, name=sequence_id, description=f"Mock {sequence_id} vector")
            
            # Add required annotations
            record.annotations = {'molecule_type': 'DNA', 'topology': 'circular'}
            
            # Add some features
            ori_feature = SeqFeature(
                FeatureLocation(100, 700, 1),
                type="rep_origin",
                qualifiers={"label": ["ori"], "note": ["Origin of replication"]}
            )
            
            amp_feature = SeqFeature(
                FeatureLocation(1000, 1861, 1),
                type="CDS",
                qualifiers={"label": ["ampR"], "gene": ["ampR"], "product": ["beta-lactamase"]}
            )
            
            mcs_feature = SeqFeature(
                FeatureLocation(2000, 2100, 1),
                type="misc_feature",
                qualifiers={"label": ["MCS"], "note": ["Multiple Cloning Site"]}
            )
            
            record.features = [ori_feature, amp_feature, mcs_feature]
            
        elif sequence_type == "gene":
            # Create a mock gene sequence (about 1kb)
            seq = Seq("ATGGTAAGCAAGGGCGAGGAG" + "NNNNNNNNNNNNNNNNNNNN" * 45 + "TAATGA")  # ~1kb
            record = SeqRecord(seq, id=sequence_id, name=sequence_id, description=f"Mock {sequence_id} gene")
            
            # Add required annotations
            record.annotations = {'molecule_type': 'DNA'}
            
            # Add coding sequence feature
            cds_feature = SeqFeature(
                FeatureLocation(0, len(seq), 1),
                type="CDS",
                qualifiers={"label": [sequence_id], "gene": [sequence_id]}
            )
            
            record.features = [cds_feature]
            
        else:
            # Generic sequence
            seq = Seq("ATGCATGCATGCATGCATGCATGCATGCATGCATGC" * 10)
            record = SeqRecord(seq, id=sequence_id, name=sequence_id, description=f"Mock {sequence_id} sequence")
            record.annotations = {'molecule_type': 'DNA'}
        
        return record
    
    def _add_promoter(self, record: SeqRecord, promoter_id: str, insertion_site: int) -> None:
        """
        Add a promoter feature to the construct
        
        Args:
            record: SeqRecord to modify
            promoter_id: ID of the promoter
            insertion_site: Position of the insert
        """
        # Try to load the actual promoter sequence
        promoter = get_sequence(promoter_id)
        
        # Determine promoter length and position
        if promoter:
            promoter_length = len(promoter.seq)
            promoter_start = max(0, insertion_site - promoter_length)
        else:
            # Use a default length if not found
            promoter_length = 100
            promoter_start = max(0, insertion_site - promoter_length)
        
        # Create the promoter feature
        promoter_feature = SeqFeature(
            FeatureLocation(promoter_start, insertion_site, 1),
            type="promoter",
            qualifiers={"label": [promoter_id], "gene": [promoter_id]}
        )
        
        # Add to features
        record.features.append(promoter_feature)
    
    def _add_feature(self, record: SeqRecord, feature_info: Dict[str, Any]) -> None:
        """
        Add a custom feature to the construct
        
        Args:
            record: SeqRecord to modify
            feature_info: Dictionary with feature information
        """
        feature_name = feature_info.get("name", "unknown_feature")
        feature_type = feature_info.get("type", "misc_feature")
        
        # For a real implementation, we'd determine the actual position
        # Here we'll just add it after the last feature or at position 100
        if record.features:
            last_feature = max(record.features, key=lambda f: f.location.end)
            start_pos = last_feature.location.end + 10
        else:
            start_pos = 100
        
        # Create a reasonable length based on feature type
        if feature_type == "tag":
            length = 30  # Tags are usually short
        elif feature_type == "promoter":
            length = 100
        elif feature_type == "terminator":
            length = 50
        else:
            length = 20
        
        # Create the feature
        new_feature = SeqFeature(
            FeatureLocation(start_pos, start_pos + length, 1),
            type=feature_type,
            qualifiers={"label": [feature_name], "note": [f"Added {feature_type}: {feature_name}"]}
        )
        
        # Add to features
        record.features.append(new_feature)

def create_construct(construct_info):
    """
    Standalone function to create a construct from provided information.
    This serves as an adapter for the app.py file to use.
    
    Args:
        construct_info: Dictionary with construct specifications
        
    Returns:
        Dictionary with construct results
    """
    # Create an instance of ConstructBuilder
    builder = ConstructBuilder()
    
    # Call the create_construct method on the instance
    return builder.create_construct(construct_info) 