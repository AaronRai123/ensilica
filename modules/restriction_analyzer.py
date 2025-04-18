"""
Restriction site analyzer for DNA sequences
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import *
from typing import Dict, Any, List, Tuple
import pandas as pd
import io
import matplotlib.pyplot as plt
import base64
import logging

def analyze_restriction_sites(record: SeqRecord, enzymes: List[str] = None) -> Dict[str, Any]:
    """
    Analyze restriction sites in a DNA sequence
    
    Args:
        record: SeqRecord to analyze
        enzymes: List of enzyme names to analyze (optional)
        
    Returns:
        Dictionary containing analysis results
    """
    # Create a RestrictionBatch with the specified enzymes
    try:
        if enzymes:
            # Filter out invalid enzymes and provide warning
            valid_enzymes = []
            invalid_enzymes = []
            for e in enzymes:
                if hasattr(Restriction, e):
                    valid_enzymes.append(getattr(Restriction, e))
                else:
                    invalid_enzymes.append(e)
            
            if invalid_enzymes:
                logger = logging.getLogger(__name__)
                logger.warning(f"Invalid enzyme names ignored: {', '.join(invalid_enzymes)}")
            
            if not valid_enzymes:
                # No valid enzymes, return empty results
                return {
                    "sites": [],
                    "dataframe": pd.DataFrame(),
                    "unique_cutters": [],
                    "rare_cutters": [],
                    "map_image": None
                }
            
            rb = RestrictionBatch(valid_enzymes)
        else:
            # Use common restriction enzymes if none specified
            common_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NdeI', 'XbaI', 
                             'PstI', 'SalI', 'SmaI', 'KpnI', 'SacI', 'BglII']
            valid_enzymes = []
            for e in common_enzymes:
                if hasattr(Restriction, e):
                    valid_enzymes.append(getattr(Restriction, e))
            
            rb = RestrictionBatch(valid_enzymes)
        
        # Search for restriction sites
        result = rb.search(record.seq)
        
        # Format the results
        sites_data = []
        for enzyme, positions in result.items():
            for pos in positions:
                # Get recognition site and cut positions
                recog_site = enzyme.site
                fwd_cut = enzyme.fst3
                rev_cut = enzyme.scd3
                
                # Calculate cut position in the sequence
                cut_pos = pos + fwd_cut if fwd_cut > 0 else pos + len(recog_site) + fwd_cut
                
                # Add to results
                sites_data.append({
                    "enzyme": enzyme.__name__,
                    "position": pos + 1,  # 1-indexed position
                    "recognition_site": recog_site,
                    "cut_position": cut_pos + 1,  # 1-indexed position
                    "overhang": enzyme.ovhg
                })
        
        # Create dataframe for easy analysis
        sites_df = pd.DataFrame(sites_data)
        
        # Sort by position
        if not sites_df.empty:
            sites_df = sites_df.sort_values("position")
        
        # Analyze for useful restriction sites
        unique_cutters = []
        rare_cutters = []
        
        if not sites_df.empty:
            enzyme_counts = sites_df["enzyme"].value_counts()
            
            # Find unique cutters (cut only once)
            unique_cutters = enzyme_counts[enzyme_counts == 1].index.tolist()
            
            # Find rare cutters (cut 2-3 times)
            rare_cutters = enzyme_counts[(enzyme_counts >= 2) & (enzyme_counts <= 3)].index.tolist()
        
        # Generate a restriction map visualization
        map_image = None
        if not sites_df.empty:
            map_image = _generate_restriction_map(record, sites_df)
        
        # Return results
        return {
            "sites": sites_data,
            "dataframe": sites_df,
            "unique_cutters": unique_cutters,
            "rare_cutters": rare_cutters,
            "map_image": map_image
        }
    except Exception as e:
        logger = logging.getLogger(__name__)
        logger.error(f"Error analyzing restriction sites: {e}")
        return {
            "sites": [],
            "dataframe": pd.DataFrame(),
            "unique_cutters": [],
            "rare_cutters": [],
            "map_image": None
        }

def _generate_restriction_map(record: SeqRecord, sites_df: pd.DataFrame) -> str:
    """
    Generate a restriction map visualization
    
    Args:
        record: SeqRecord to visualize
        sites_df: DataFrame with restriction sites
        
    Returns:
        Base64-encoded PNG image
    """
    # Create figure and axis
    fig, ax = plt.figure(figsize=(10, 5), dpi=100), plt.gca()
    
    # Draw the DNA sequence as a line
    seq_length = len(record.seq)
    ax.plot([0, seq_length], [0, 0], 'k-', linewidth=2)
    
    # Add markers for restriction sites
    positions = sites_df["position"].tolist()
    enzymes = sites_df["enzyme"].tolist()
    
    # Set colors for different enzymes
    enzyme_colors = {}
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    for i, enzyme in enumerate(set(enzymes)):
        enzyme_colors[enzyme] = colors[i % len(colors)]
    
    # Draw markers and labels
    for pos, enzyme in zip(positions, enzymes):
        # Draw marker
        ax.plot([pos, pos], [-0.1, 0.1], '-', color=enzyme_colors[enzyme], linewidth=1.5)
        
        # Draw label
        ax.text(pos, 0.15 if pos % 2 == 0 else -0.15, enzyme, 
                ha='center', va='center', rotation=90 if pos % 2 == 0 else 270,
                color=enzyme_colors[enzyme], fontsize=8)
    
    # Add scale bar
    scale_interval = seq_length // 5
    for i in range(0, seq_length + 1, scale_interval):
        ax.plot([i, i], [-0.05, 0.05], 'k-', linewidth=1)
        ax.text(i, -0.1, f"{i}", ha='center', va='center', fontsize=8)
    
    # Add title and labels
    ax.set_title(f"Restriction Map for {record.id}")
    ax.set_xlabel("Base Position")
    ax.set_xlim(-100, seq_length + 100)
    ax.set_ylim(-0.5, 0.5)
    ax.axis('off')
    
    # Save figure to a buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close(fig)
    
    # Encode image as base64
    buf.seek(0)
    map_image = base64.b64encode(buf.read()).decode('utf-8')
    
    return map_image

def get_compatible_enzyme_pairs(record: SeqRecord) -> List[Dict[str, Any]]:
    """
    Find compatible enzyme pairs for cloning
    
    Args:
        record: SeqRecord to analyze
        
    Returns:
        List of dictionaries with compatible enzyme pairs
    """
    try:
        # Get common restriction enzymes
        common_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NdeI', 'XbaI', 
                          'PstI', 'SalI', 'SmaI', 'KpnI', 'SacI', 'BglII',
                          'NotI', 'XmaI', 'SpeI', 'SphI', 'EcoRV', 'NcoI',
                          'ApaI', 'BstXI', 'ClaI', 'AgeI', 'SacII', 'SnaBI']
        
        # Create RestrictionBatch with valid enzymes only
        valid_enzymes = []
        for e in common_enzymes:
            if hasattr(Restriction, e):
                valid_enzymes.append(getattr(Restriction, e))
        
        # If no valid enzymes, return empty list
        if not valid_enzymes:
            return []
            
        rb = RestrictionBatch(valid_enzymes)
        
        # Analyze the sequence
        result = rb.search(record.seq)
        
        # Find unique cutters
        unique_cutters = []
        for enzyme, positions in result.items():
            if len(positions) == 1:
                unique_cutters.append((enzyme, positions[0]))
        
        # Find compatible pairs
        pairs = []
        for i in range(len(unique_cutters)):
            for j in range(i+1, len(unique_cutters)):
                enzyme1, pos1 = unique_cutters[i]
                enzyme2, pos2 = unique_cutters[j]
                
                # Ensure enzymes cut at different positions
                if pos1 != pos2:
                    # Ensure proper order for fragment isolation
                    if pos1 > pos2:
                        enzyme1, pos1, enzyme2, pos2 = enzyme2, pos2, enzyme1, pos1
                    
                    # Calculate fragment size
                    fragment_size = pos2 - pos1
                    
                    # Check if the fragment is a reasonable size (50bp to 1/2 of sequence)
                    if 50 <= fragment_size <= len(record.seq) // 2:
                        pairs.append({
                            "enzyme1": enzyme1.__name__,
                            "enzyme2": enzyme2.__name__,
                            "position1": pos1 + 1,  # 1-indexed
                            "position2": pos2 + 1,  # 1-indexed
                            "fragment_size": fragment_size,
                            "compatibility": _check_buffer_compatibility(enzyme1, enzyme2)
                        })
        
        # Sort by compatibility and fragment size
        pairs.sort(key=lambda x: (-x["compatibility"], x["fragment_size"]))
        
        return pairs
    except Exception as e:
        logger = logging.getLogger(__name__)
        logger.error(f"Error finding compatible enzyme pairs: {e}")
        return []

def _check_buffer_compatibility(enzyme1, enzyme2) -> int:
    """
    Check buffer compatibility between two enzymes
    
    Args:
        enzyme1: First restriction enzyme
        enzyme2: Second restriction enzyme
        
    Returns:
        Compatibility score (0-10)
    """
    # This is a simplified model of compatibility
    # In reality, you would need to check actual buffer compositions
    
    # Optimal temperature difference - closer is better
    temp_diff = abs(enzyme1.opt_temp - enzyme2.opt_temp)
    
    # Optimal buffer differences
    buffer_diff = 5  # Default medium compatibility
    
    # Some known compatible pairs
    compatible_pairs = [
        ("EcoRI", "HindIII"),
        ("BamHI", "HindIII"),
        ("EcoRI", "BamHI"),
        ("XhoI", "SalI"),
        ("NdeI", "XhoI"),
        ("XbaI", "SpeI"),
        ("BglII", "BamHI")
    ]
    
    # Check if this is a known compatible pair
    if (enzyme1.__name__, enzyme2.__name__) in compatible_pairs or \
       (enzyme2.__name__, enzyme1.__name__) in compatible_pairs:
        buffer_diff = 0  # Perfect compatibility
    
    # Calculate overall compatibility score (0-10 scale)
    # Lower temp difference and buffer difference = higher score
    score = 10 - min(5, temp_diff // 5) - min(5, buffer_diff)
    
    return score 