"""
Primer Tools module for designing primers using Primer3
"""
from typing import Dict, Any, List
import os
import logging

# Initialize logger
logger = logging.getLogger("primer_tools")
logger.setLevel(logging.INFO)

try:
    import primer3
    PRIMER3_AVAILABLE = True
except ImportError:
    PRIMER3_AVAILABLE = False
    print("Warning: primer3-py not installed. Using simplified calculations.")

def design_primers(insert_region: Dict[str, Any], construct_info: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design primers for the inserted gene using simple method since Primer3 is failing.
    
    Args:
        insert_region: Dictionary containing information about the inserted region
        construct_info: Dictionary containing information about the construct
        
    Returns:
        Dictionary containing primer information
        
    Raises:
        ValueError: If primer design fails or sequence is invalid
    """
    gene_name = construct_info.get("gene_name", "unknown_gene")
    template_seq = insert_region["sequence"]
    cloning_method = construct_info.get("cloning_method", "Gibson Assembly")
    
    # Validate the sequence
    if not template_seq:
        raise ValueError(f"Could not generate primers: empty sequence for {gene_name}")
    
    if len(template_seq) < 60:
        raise ValueError(f"Could not generate primers: sequence for {gene_name} is too short (minimum 60 bp required)")
    
    # Count ambiguous bases (N)
    seq_str = str(template_seq).upper()
    n_count = seq_str.count('N')
    n_percentage = (n_count / len(seq_str)) * 100
    
    if n_percentage > 10:
        raise ValueError(f"Could not generate primers: sequence for {gene_name} contains {n_percentage:.1f}% ambiguous bases (maximum 10% allowed)")
    
    # Determine if we need to add restriction sites
    restriction_sites = None
    if cloning_method == "Restriction Enzyme Cloning":
        restriction_sites = get_restriction_sites_for_cloning(construct_info)
    
    try:
        # Use simple primer design instead of Primer3
        logger.info(f"Designing simple primers for {gene_name} (length: {len(template_seq)})")
        primers = design_simple_primers(template_seq, gene_name, restriction_sites)
        
        return {
            "primers": primers,
            "target_region": {
                "start": insert_region["start"],
                "end": insert_region["end"],
                "length": insert_region["end"] - insert_region["start"]
            }
        }
    except Exception as e:
        error_msg = f"Could not generate primers for {gene_name}: {str(e)}"
        logger.error(error_msg)
        raise ValueError(error_msg) from e

def get_restriction_sites_for_cloning(construct_info: Dict[str, Any]) -> Dict[str, Dict[str, str]]:
    """
    Get appropriate restriction sites for cloning based on construct info
    
    Args:
        construct_info: Dictionary containing construct information
        
    Returns:
        Dictionary with 'forward' and 'reverse' keys for restriction sites
    """
    # Common restriction sites with overhangs
    sites = {
        "EcoRI": "GAATTC",
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
        "XhoI": "CTCGAG",
        "NdeI": "CATATG",
        "XbaI": "TCTAGA"
    }
    
    # Default sites
    default_sites = {
        "forward": {"name": "EcoRI", "sequence": "GAATTC", "overhang": "G"},
        "reverse": {"name": "HindIII", "sequence": "AAGCTT", "overhang": "G"}
    }
    
    # Check if construct_info specifies restriction enzymes
    if "restriction_enzymes" in construct_info and isinstance(construct_info["restriction_enzymes"], list):
        enzymes = construct_info["restriction_enzymes"]
        if len(enzymes) >= 2:
            # Use the specified enzymes
            if enzymes[0] in sites and enzymes[1] in sites:
                return {
                    "forward": {"name": enzymes[0], "sequence": sites[enzymes[0]], "overhang": "G"},
                    "reverse": {"name": enzymes[1], "sequence": sites[enzymes[1]], "overhang": "G"}
                }
    
    return default_sites

def design_primers_with_primer3(sequence: str, flank_size: int, template_length: int, gene_name: str, 
                               restriction_sites: Dict[str, Dict[str, str]] = None) -> List[Dict[str, Any]]:
    """
    Design primers using Primer3 with optional restriction sites
    
    Args:
        sequence: Full sequence including flanking regions
        flank_size: Size of flanking regions
        template_length: Length of the actual template (without flanks)
        gene_name: Name of the gene for primer names
        restriction_sites: Optional restriction sites to add
        
    Returns:
        List of dictionaries with primer information
    """
    # Define the target region (excluding flanks)
    target_start = flank_size
    target_length = template_length
    
    # Setup Primer3 parameters
    seq_args = {
        'SEQUENCE_ID': gene_name,
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': [0, len(sequence)],
        'SEQUENCE_TARGET': [target_start, target_length],
    }
    
    # Primer design parameters - make more permissive
    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 30,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_MIN_GC': 30.0,
        'PRIMER_MAX_GC': 70.0,
        'PRIMER_NUM_RETURN': 3,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_MAX_END_STABILITY': 9.0,
        'PRIMER_EXPLAIN_FLAG': 1,
        'PRIMER_MAX_POLY_X': 5,
        'PRIMER_LIBERAL_BASE': 1
    }
    
    # Run Primer3
    primer3_results = primer3.bindings.designPrimers(seq_args, global_args)
    
    # Extract primer information
    primers = []
    
    # Forward primer
    if 'PRIMER_LEFT_0_SEQUENCE' in primer3_results:
        forward_seq = primer3_results['PRIMER_LEFT_0_SEQUENCE']
        forward_tm = primer3_results['PRIMER_LEFT_0_TM']
        forward_gc = primer3_results['PRIMER_LEFT_0_GC_PERCENT']
        
        # Add restriction site if needed
        if restriction_sites:
            fwd_site = restriction_sites['forward']
            forward_seq_with_site = fwd_site['overhang'] + fwd_site['sequence'] + forward_seq
            forward_props = calculate_primer_properties(forward_seq_with_site)
            
            forward_primer = {
                "name": f"{gene_name}_F",
                "sequence": forward_seq_with_site,
                "tm": round(forward_props["tm"], 1),
                "gc_content": round(forward_props["gc_content"], 1),
                "length": len(forward_seq_with_site),
                "restriction_site": fwd_site['name']
            }
        else:
            forward_primer = {
                "name": f"{gene_name}_F",
                "sequence": forward_seq,
                "tm": round(forward_tm, 1),
                "gc_content": round(forward_gc, 1),
                "length": len(forward_seq)
            }
        
        primers.append(forward_primer)
    
    # Reverse primer
    if 'PRIMER_RIGHT_0_SEQUENCE' in primer3_results:
        reverse_seq = primer3_results['PRIMER_RIGHT_0_SEQUENCE']
        reverse_tm = primer3_results['PRIMER_RIGHT_0_TM']
        reverse_gc = primer3_results['PRIMER_RIGHT_0_GC_PERCENT']
        
        # Add restriction site if needed
        if restriction_sites:
            rev_site = restriction_sites['reverse']
            reverse_seq_with_site = rev_site['overhang'] + rev_site['sequence'] + reverse_seq
            reverse_props = calculate_primer_properties(reverse_seq_with_site)
            
            reverse_primer = {
                "name": f"{gene_name}_R",
                "sequence": reverse_seq_with_site,
                "tm": round(reverse_props["tm"], 1),
                "gc_content": round(reverse_props["gc_content"], 1),
                "length": len(reverse_seq_with_site),
                "restriction_site": rev_site['name']
            }
        else:
            reverse_primer = {
                "name": f"{gene_name}_R",
                "sequence": reverse_seq,
                "tm": round(reverse_tm, 1),
                "gc_content": round(reverse_gc, 1),
                "length": len(reverse_seq)
            }
        
        primers.append(reverse_primer)
    
    # If Primer3 failed to generate both primers, raise an error
    if len(primers) < 2:
        raise ValueError(f"Primer3 failed to design suitable primers for {gene_name}")
    
    # For verification primers, design internal primers
    template = sequence[flank_size:flank_size+template_length]
    verify_primers = design_verification_primers(template, gene_name)
    
    # Add verification primers to the list
    primers.extend(verify_primers)
    
    return primers

def design_fixed_primers(template: str, gene_name: str) -> List[Dict[str, Any]]:
    """Design primers at fixed positions (start and end of template)."""
    # Forward primer: first 20-25 bp
    forward_length = min(25, len(template) // 3)
    forward_seq = template[:forward_length]
    
    # Reverse primer: last 20-25 bp (reverse complement)
    reverse_length = min(25, len(template) // 3)
    reverse_template = template[-reverse_length:]
    reverse_seq = reverse_complement(reverse_template)
    
    # Calculate properties
    forward_props = calculate_primer_properties(forward_seq)
    reverse_props = calculate_primer_properties(reverse_seq)
    
    forward_primer = {
        "name": f"{gene_name}_F",
        "sequence": forward_seq,
        "tm": forward_props["tm"],
        "gc_content": forward_props["gc_content"],
        "length": len(forward_seq)
    }
    
    reverse_primer = {
        "name": f"{gene_name}_R",
        "sequence": reverse_seq,
        "tm": reverse_props["tm"],
        "gc_content": reverse_props["gc_content"],
        "length": len(reverse_seq)
    }
    
    # Create verification primers as well (slightly inset)
    verify_offset = min(30, len(template) // 4)
    
    verification_f_seq = template[verify_offset:verify_offset+25]
    verification_f_props = calculate_primer_properties(verification_f_seq)
    
    verification_r_template = template[-(verify_offset+25):-verify_offset]
    verification_r_seq = reverse_complement(verification_r_template)
    verification_r_props = calculate_primer_properties(verification_r_seq)
    
    verification_f = {
        "name": f"{gene_name}_verify_F",
        "sequence": verification_f_seq,
        "tm": verification_f_props["tm"],
        "gc_content": verification_f_props["gc_content"],
        "length": len(verification_f_seq)
    }
    
    verification_r = {
        "name": f"{gene_name}_verify_R",
        "sequence": verification_r_seq,
        "tm": verification_r_props["tm"],
        "gc_content": verification_r_props["gc_content"],
        "length": len(verification_r_seq)
    }
    
    return [forward_primer, reverse_primer, verification_f, verification_r]

def design_simple_primers(template: str, gene_name: str, restriction_sites: Dict[str, Dict[str, str]] = None) -> List[Dict[str, Any]]:
    """Design primers using a simple approach (for fallback)."""
    # Forward primer: first 20-25 bp
    forward_length = min(25, max(18, len(template) // 10))
    forward_seq = template[:forward_length]
    
    # Reverse primer: last 20-25 bp (reverse complement)
    reverse_length = min(25, max(18, len(template) // 10))
    reverse_template = template[-reverse_length:]
    reverse_seq = reverse_complement(reverse_template)
    
    # Calculate properties
    forward_props = calculate_primer_properties(forward_seq)
    reverse_props = calculate_primer_properties(reverse_seq)
    
    forward_primer = {
        "name": f"{gene_name}_F",
        "sequence": forward_seq,
        "tm": forward_props["tm"],
        "gc_content": forward_props["gc_content"],
        "length": len(forward_seq)
    }
    
    reverse_primer = {
        "name": f"{gene_name}_R",
        "sequence": reverse_seq,
        "tm": reverse_props["tm"],
        "gc_content": reverse_props["gc_content"],
        "length": len(reverse_seq)
    }
    
    # Add verification primers
    verify_primers = design_verification_primers(template, gene_name)
    
    return [forward_primer, reverse_primer] + verify_primers

def design_verification_primers(template: str, gene_name: str) -> List[Dict[str, Any]]:
    """Design internal verification primers."""
    # Only design verification primers if the template is long enough
    if len(template) < 100:
        return []
    
    # Forward verification primer (start 30bp in)
    offset = min(30, len(template) // 5)
    verify_f_length = min(25, max(18, len(template) // 15))
    verify_f_seq = template[offset:offset+verify_f_length]
    
    # Reverse verification primer (end 30bp from the end)
    verify_r_template = template[-(offset+verify_f_length):-offset]
    verify_r_seq = reverse_complement(verify_r_template)
    
    # Calculate properties
    verify_f_props = calculate_primer_properties(verify_f_seq)
    verify_r_props = calculate_primer_properties(verify_r_seq)
    
    verification_f = {
        "name": f"{gene_name}_verify_F",
        "sequence": verify_f_seq,
        "tm": verify_f_props["tm"],
        "gc_content": verify_f_props["gc_content"],
        "length": len(verify_f_seq)
    }
    
    verification_r = {
        "name": f"{gene_name}_verify_R",
        "sequence": verify_r_seq,
        "tm": verify_r_props["tm"],
        "gc_content": verify_r_props["gc_content"],
        "length": len(verify_r_seq)
    }
    
    return [verification_f, verification_r]

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(seq))

def calculate_primer_properties(sequence: str) -> Dict[str, Any]:
    """
    Calculate properties of a primer sequence using Biopython
    
    Args:
        sequence: Primer sequence (5' to 3')
        
    Returns:
        Dictionary with Tm, GC content, and other properties
    """
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq
    
    seq = Seq(sequence)
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    
    # Calculate melting temperature using nearest-neighbor method
    tm = mt.Tm_NN(seq, nn_table=mt.DNA_NN4, Na=50, Mg=1.5, dNTPs=0.6)
    
    return {
        "tm": tm,
        "gc_content": gc_content,
        "length": len(sequence)
    } 