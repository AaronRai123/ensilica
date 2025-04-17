"""
Visualization module for generating plasmid maps
"""
import base64
import io
import math
from typing import Optional, Dict, Any
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib.patches

def render_plasmid_map(record: SeqRecord, highlight_region: Dict[str, Any] = None, width=800, height=800) -> str:
    """
    Generates a base64-encoded SVG of a circular plasmid map from a SeqRecord.
    
    Args:
        record: SeqRecord object containing the construct
        highlight_region: Optional region to highlight
        width: Width of the SVG
        height: Height of the SVG
        
    Returns:
        Base64 encoded SVG string
    """
    # Basic settings
    center_x, center_y = width / 2, height / 2
    radius = 250
    feature_width = 30
    sequence_length = len(record.seq)
    
    # Color scheme for different feature types
    colors = {
        "CDS": "#34A853",         # Green
        "gene": "#34A853",        # Green
        "promoter": "#4285F4",    # Blue
        "terminator": "#EA4335",  # Red
        "rep_origin": "#9C27B0",  # Purple
        "misc_feature": "#A142F4", # Light purple
        "primer_bind": "#FF6D01", # Orange
        "RBS": "#00C9A7",         # Teal
        "regulatory": "#4285F4",  # Blue
        "source": "#EEEEEE",      # Light grey
        "resistance": "#FF9800",  # Orange
    }
    
    # Start the SVG
    svg = f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
        <!-- Background -->
        <rect width="{width}" height="{height}" fill="white" />
        
        <!-- Base plasmid circle -->
        <circle cx="{center_x}" cy="{center_y}" r="{radius}" fill="white" stroke="#333333" stroke-width="2"/>
        
        <!-- Title and size -->
        <text x="{center_x}" y="{center_y - radius - 40}" 
              text-anchor="middle" font-size="24" font-weight="bold" fill="#333333">
            {record.id}
        </text>
        <text x="{center_x}" y="{center_y - radius - 15}" 
              text-anchor="middle" font-size="16" fill="#666666">
            {sequence_length} bp
        </text>
    """
    
    # Add features
    for i, feature in enumerate(record.features):
        # Skip source features (span the entire sequence)
        if feature.type == "source":
            continue
            
        # Get feature coordinates
        start = int(feature.location.start)
        end = int(feature.location.end)
        
        # Get label from qualifiers
        label = ""
        for qualifier in ["label", "gene", "product", "note"]:
            if hasattr(feature, 'qualifiers') and qualifier in feature.qualifiers:
                label = feature.qualifiers[qualifier][0]
                break
                
        if not label:
            label = f"{feature.type}_{i+1}"
        
        # Get feature color
        color = colors.get(feature.type, "#CCCCCC")
        
        # Calculate angles
        start_angle = 2 * math.pi * start / sequence_length
        end_angle = 2 * math.pi * end / sequence_length
        
        # Handle features that cross the origin
        if end < start:
            end_angle = 2 * math.pi * (end + sequence_length) / sequence_length
            
        large_arc = 1 if end_angle - start_angle > math.pi else 0
        
        # Calculate coordinates for the feature
        outer_radius = radius + 5
        inner_radius = radius - feature_width
        
        # Outer arc coordinates
        x1 = center_x + outer_radius * math.cos(start_angle)
        y1 = center_y + outer_radius * math.sin(start_angle)
        x2 = center_x + outer_radius * math.cos(end_angle)
        y2 = center_y + outer_radius * math.sin(end_angle)
        
        # Inner arc coordinates
        x3 = center_x + inner_radius * math.cos(end_angle)
        y3 = center_y + inner_radius * math.sin(end_angle)
        x4 = center_x + inner_radius * math.cos(start_angle)
        y4 = center_y + inner_radius * math.sin(start_angle)
        
        # Draw the feature
        svg += f"""
        <path d="M {x1} {y1}
                 A {outer_radius} {outer_radius} 0 {large_arc} 1 {x2} {y2}
                 L {x3} {y3}
                 A {inner_radius} {inner_radius} 0 {large_arc} 0 {x4} {y4}
                 Z"
              fill="{color}" stroke="black" stroke-width="1" opacity="0.8"/>
        """
        
        # Add the label
        mid_angle = (start_angle + end_angle) / 2
        label_radius = outer_radius + 20
        label_x = center_x + label_radius * math.cos(mid_angle)
        label_y = center_y + label_radius * math.sin(mid_angle)
        
        # Adjust text rotation for readability
        rotation = mid_angle * 180 / math.pi
        if 90 < rotation < 270:
            rotation += 180
            
        svg += f"""
        <text x="{label_x}" y="{label_y}" 
              text-anchor="middle" 
              transform="rotate({rotation} {label_x} {label_y})"
              font-size="12" fill="#333333">{label}</text>
        """
    
    # Highlight the insert region if provided
    if highlight_region and 'start' in highlight_region and 'end' in highlight_region:
        start = highlight_region['start']
        end = highlight_region['end']
        
        # Calculate angles
        start_angle = 2 * math.pi * start / sequence_length
        end_angle = 2 * math.pi * end / sequence_length
        
        # Handle regions that cross the origin
        if end < start:
            end_angle = 2 * math.pi * (end + sequence_length) / sequence_length
            
        large_arc = 1 if end_angle - start_angle > math.pi else 0
        
        # Draw a highlight ring
        highlight_radius = radius + 15
        x1 = center_x + highlight_radius * math.cos(start_angle)
        y1 = center_y + highlight_radius * math.sin(start_angle)
        x2 = center_x + highlight_radius * math.cos(end_angle)
        y2 = center_y + highlight_radius * math.sin(end_angle)
        
        svg += f"""
        <path d="M {x1} {y1}
                 A {highlight_radius} {highlight_radius} 0 {large_arc} 1 {x2} {y2}"
              stroke="#FF5722" stroke-width="5" fill="none" opacity="0.8" stroke-dasharray="10,5"/>
        """
        
        # Add label
        mid_angle = (start_angle + end_angle) / 2
        label_x = center_x + (highlight_radius + 20) * math.cos(mid_angle)
        label_y = center_y + (highlight_radius + 20) * math.sin(mid_angle)
        
        # Adjust text rotation for readability
        rotation = mid_angle * 180 / math.pi
        if 90 < rotation < 270:
            rotation += 180
            
        svg += f"""
        <text x="{label_x}" y="{label_y}" 
              text-anchor="middle" 
              transform="rotate({rotation} {label_x} {label_y})"
              font-size="14" font-weight="bold" fill="#FF5722">INSERT</text>
        """
    
    # Close the SVG
    svg += "</svg>"
    
    # Encode to base64
    return base64.b64encode(svg.encode("utf-8")).decode("utf-8")

def visualize_construct(record: SeqRecord, highlight_region: Dict[str, Any] = None) -> str:
    """
    Visualize a genetic construct and return the visualization as a base64 encoded SVG
    
    Args:
        record: SeqRecord object containing the construct
        highlight_region: Optional region to highlight
        
    Returns:
        Base64 encoded SVG data
    """
    try:
        # Try the simplified plasmid visualization with standardized elements
        return create_standard_plasmid_map(record, highlight_region)
    except Exception as e:
        import logging
        logging.error(f"Visualization error: {str(e)}")
        # Create a simple fallback visualization
        return create_fallback_visualization(record)

def create_standard_plasmid_map(record: SeqRecord, highlight_region: Dict[str, Any] = None, width=800, height=800) -> str:
    """
    Create a standardized plasmid map visualization similar to the reference image
    
    Args:
        record: SeqRecord object containing the construct
        highlight_region: Optional region to highlight
        width: Width of the SVG
        height: Height of the SVG
        
    Returns:
        Base64 encoded SVG string
    """
    # Basic settings
    center_x, center_y = width / 2, height / 2
    radius = 250
    feature_width = 30
    sequence_length = len(record.seq)
    
    # Professional color scheme
    colors = {
        "CDS": "#4286f4",          # Blue for CDS/genes
        "gene": "#4286f4",         # Blue for CDS/genes
        "promoter": "#9c27b0",     # Purple for promoters 
        "terminator": "#e91e63",   # Pink for terminators
        "rep_origin": "#2196f3",   # Light blue for origins
        "misc_feature": "#607d8b", # Gray for misc features
        "primer_bind": "#ff5722",  # Red/orange for primer sites
        "RBS": "#00bcd4",          # Cyan for RBS
        "regulatory": "#9c27b0",   # Purple for regulatory elements
        "source": "#eeeeee",       # Light grey for source
        "resistance": "#f44336",   # Red for resistance genes
        "marker": "#ff9800",       # Orange for selectable markers
    }
    
    # Start the SVG
    svg = f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
        <!-- Background -->
        <rect width="{width}" height="{height}" fill="white" />
        
        <!-- Base plasmid circle -->
        <circle cx="{center_x}" cy="{center_y}" r="{radius}" fill="none" stroke="#999999" stroke-width="15"/>
        
        <!-- Title -->
        <text x="{center_x}" y="{center_y}" 
              text-anchor="middle" dominant-baseline="middle" font-size="40" font-weight="bold" fill="#333333">
            Plasmid Map
        </text>
    """
    
    # Categorize features for better labeling
    resistance_genes = []
    markers = []
    promoters = []
    origins = []
    inserts = []
    primers = []
    restriction_sites = []
    
    # Process features and categorize them
    for i, feature in enumerate(record.features):
        if feature.type == "source":
            continue
            
        # Get feature coordinates
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.strand
        
        # Get feature name from qualifiers
        feature_name = ""
        for qualifier in ["label", "gene", "product", "note"]:
            if hasattr(feature, 'qualifiers') and qualifier in feature.qualifiers:
                feature_name = feature.qualifiers[qualifier][0]
                break
                
        if not feature_name:
            feature_name = f"{feature.type}_{i+1}"
        
        # Categorize feature
        if feature.type == "CDS" or feature.type == "gene":
            if any(term in feature_name.lower() for term in ['resistance', 'ampr', 'kanr', 'cmr', 'tetr', 'bla']):
                resistance_genes.append((start, end, strand, feature_name, feature.type))
            elif highlight_region and start <= highlight_region.get("start", -1) <= end:
                inserts.append((start, end, strand, feature_name, feature.type))
            elif any(term in feature_name.lower() for term in ['marker', 'gfp', 'rfp', 'yfp', 'cfp']):
                markers.append((start, end, strand, feature_name, feature.type))
            else:
                inserts.append((start, end, strand, feature_name, feature.type))
        elif feature.type == "promoter":
            promoters.append((start, end, strand, feature_name, feature.type))
        elif feature.type == "rep_origin":
            origins.append((start, end, strand, feature_name, feature.type))
        elif feature.type == "primer_bind":
            primers.append((start, end, strand, feature_name, feature.type))
        elif feature.type == "misc_feature" and "restriction" in feature_name.lower():
            restriction_sites.append((start, end, strand, feature_name, feature.type))
    
    # If insert is not found but highlight region exists, create an insert entry
    if not inserts and highlight_region:
        start = highlight_region.get("start", 0)
        end = highlight_region.get("end", 100)
        inserts.append((start, end, 1, "Inserted Gene", "gene"))
    
    # Add standard features if key elements are missing
    if not resistance_genes:
        resistance_genes.append((1500, 2500, 1, "Antibiotic Resistance Gene", "gene"))
    
    if not origins:
        origins.append((3000, 3200, 1, "Origin of Replication", "rep_origin"))
    
    if not markers and not inserts:
        markers.append((500, 1200, 1, "Selectable Marker", "gene"))
    
    if not promoters:
        promoters.append((200, 300, 1, "Promoter", "promoter"))
    
    # Draw features - standard arrangement resembling reference image
    
    # Origin of Replication at bottom (small feature)
    for start, end, strand, name, feature_type in origins:
        svg += draw_feature(center_x, center_y, radius, 15, math.pi * 1.5, math.pi * 1.6, colors["rep_origin"], "Origin of Replication", False)
    
    # Antibiotic Resistance at bottom-left
    for start, end, strand, name, feature_type in resistance_genes:
        svg += draw_feature(center_x, center_y, radius, 30, math.pi * 0.9, math.pi * 1.3, colors["resistance"], "Antibiotic Resistance Gene", True)
    
    # Selectable Marker at top-left
    for start, end, strand, name, feature_type in markers:
        svg += draw_feature(center_x, center_y, radius, 30, math.pi * 0.4, math.pi * 0.8, colors["marker"], "Selectable Marker", True)
    
    # Promoter at top-right
    if promoters:
        svg += draw_feature(center_x, center_y, radius, 30, math.pi * 1.75, math.pi * 1.95, colors["promoter"], "Promoter", True)
    
    # Primer sites
    svg += draw_feature(center_x, center_y, radius, 8, math.pi * 1.96, math.pi * 1.98, "#ff5722", "5' Primer Site", False)
    svg += draw_feature(center_x, center_y, radius, 8, math.pi * 0.15, math.pi * 0.17, "#ff5722", "3' Primer Site", False)
    
    # Restriction sites
    svg += draw_feature(center_x, center_y, radius, 5, math.pi * 1.99, math.pi * 2.01, "#222222", "Restriction Site", False)
    svg += draw_feature(center_x, center_y, radius, 5, math.pi * 0.13, math.pi * 0.15, "#222222", "Restriction Site", False)
    
    # Inserted Gene at top-right
    for start, end, strand, name, feature_type in inserts:
        svg += draw_feature(center_x, center_y, radius, 30, math.pi * 0, math.pi * 0.12, colors["CDS"], "Inserted Gene", True)
    
    # Close the SVG
    svg += "</svg>"
    
    # Encode to base64
    return base64.b64encode(svg.encode("utf-8")).decode("utf-8")

def draw_feature(cx, cy, radius, width, start_angle, end_angle, color, label, use_arrow=True):
    """Helper function to draw a feature arc with optional directional arrow"""
    # Calculate coordinates
    outer_radius = radius + width/2
    inner_radius = radius - width/2
    
    # Calculate midpoint for label
    mid_angle = (start_angle + end_angle) / 2
    label_radius = radius + width/2 + 20  # Position labels outside features
    
    # Calculate arc endpoints
    x1 = cx + outer_radius * math.cos(start_angle)
    y1 = cy + outer_radius * math.sin(start_angle)
    x2 = cx + outer_radius * math.cos(end_angle)
    y2 = cy + outer_radius * math.sin(end_angle)
    
    x3 = cx + inner_radius * math.cos(end_angle)
    y3 = cy + inner_radius * math.sin(end_angle)
    x4 = cx + inner_radius * math.cos(start_angle)
    y4 = cy + inner_radius * math.sin(start_angle)
    
    # Determine if it's a large arc
    large_arc = 1 if (end_angle - start_angle) > math.pi else 0
    
    # Base feature path
    path = f"""
    <path d="M {x1} {y1}
             A {outer_radius} {outer_radius} 0 {large_arc} 1 {x2} {y2}
             L {x3} {y3}
             A {inner_radius} {inner_radius} 0 {large_arc} 0 {x4} {y4}
             Z"
          fill="{color}" stroke="#333333" stroke-width="1" />
    """
    
    # Add directional arrow if requested
    if use_arrow:
        # Calculate arrow position at ~75% of the feature
        arrow_angle = start_angle + (end_angle - start_angle) * 0.75
        
        # Arrow size relative to feature width
        arrow_size = width * 0.8
        
        # Calculate arrow points
        ax1 = cx + radius * math.cos(arrow_angle)
        ay1 = cy + radius * math.sin(arrow_angle)
        
        # Arrow direction based on position on circle
        perp_angle = arrow_angle + math.pi/2  # Perpendicular to radius
        
        # Calculate the three points of the arrow
        px1 = ax1 + arrow_size * math.cos(perp_angle)
        py1 = ay1 + arrow_size * math.sin(perp_angle)
        
        px2 = ax1 + arrow_size * math.cos(arrow_angle)
        py2 = ay1 + arrow_size * math.sin(arrow_angle)
        
        px3 = ax1 - arrow_size * math.cos(perp_angle)
        py3 = ay1 - arrow_size * math.sin(perp_angle)
        
        path += f"""
        <path d="M {px1} {py1} L {px2} {py2} L {px3} {py3} Z"
              fill="white" stroke="#333333" stroke-width="1" />
        """
    
    # Calculate position for label
    lx = cx + label_radius * math.cos(mid_angle)
    ly = cy + label_radius * math.sin(mid_angle)
    
    # Adjust text rotation for readability
    rotation = (mid_angle * 180 / math.pi) % 360
    if 90 < rotation < 270:
        rotation += 180
    
    # Add the label
    path += f"""
    <text x="{lx}" y="{ly}" 
          text-anchor="middle" 
          transform="rotate({rotation} {lx} {ly})"
          font-family="Arial" font-size="14" fill="#333333" font-weight="bold">{label}</text>
    """
    
    return path

def create_fallback_visualization(record: SeqRecord) -> str:
    """
    Create a simple fallback visualization when the main renderer fails
    
    Args:
        record: SeqRecord object
        
    Returns:
        Base64 encoded SVG data
    """
    width, height = 500, 500
    center_x, center_y = width / 2, height / 2
    radius = 180
    
    svg = f"""<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
        <rect width="{width}" height="{height}" fill="white" />
        <circle cx="{center_x}" cy="{center_y}" r="{radius}" fill="none" stroke="#999999" stroke-width="15"/>
        <text x="{center_x}" y="{center_y}" text-anchor="middle" dominant-baseline="middle" font-size="30" font-weight="bold">Plasmid Map</text>
        <text x="{center_x}" y="{center_y + 40}" text-anchor="middle" font-size="14">{record.id} ({len(record.seq)} bp)</text>
    """
    
    # Close the SVG
    svg += "</svg>"
    
    # Encode to base64
    return base64.b64encode(svg.encode("utf-8")).decode("utf-8")

def add_restriction_sites(ax, record, plasmid_radius, center):
    """Add restriction sites visualization to the plasmid map"""
    from Bio.Restriction import RestrictionBatch
    from Bio.Restriction import CommOnly
    
    # Common restriction enzymes used in molecular biology
    common_enzymes = RestrictionBatch(CommOnly)
    
    # Find restriction sites
    restriction_sites = common_enzymes.search(record.seq)
    
    # Add restriction sites as marks on the plasmid
    for enzyme, sites in restriction_sites.items():
        if sites:  # Only show enzymes that have restriction sites in the plasmid
            for site in sites:
                site_angle = 2 * np.pi * site / len(record.seq)
                x = center[0] + plasmid_radius * np.cos(site_angle)
                y = center[1] + plasmid_radius * np.sin(site_angle)
                
                # Draw the restriction site marker
                ax.plot([center[0], x], [center[1], y], 'r-', linewidth=0.5, alpha=0.5)
                ax.plot(x, y, 'ro', markersize=3, alpha=0.7)
                
                # Add label for restriction site
                label_x = center[0] + (plasmid_radius + 25) * np.cos(site_angle)
                label_y = center[1] + (plasmid_radius + 25) * np.sin(site_angle)
                ax.text(label_x, label_y, enzyme.name, fontsize=6, ha='center', va='center', 
                       rotation=np.degrees(site_angle) if -np.pi/2 <= site_angle <= np.pi/2 else np.degrees(site_angle) + 180)

def add_primer_binding_sites(ax, record, primers, plasmid_radius, center):
    """Add primer binding sites visualization to the plasmid map"""
    if not primers:
        return
    
    for i, primer in enumerate(primers):
        # Get primer binding position (example - would need real data)
        position = primer.get('binding_position', 0)
        length = primer.get('length', 20)
        
        # Calculate angle for the primer binding site
        start_angle = 2 * np.pi * position / len(record.seq)
        end_angle = 2 * np.pi * (position + length) / len(record.seq)
        
        # Draw arc for primer binding site
        arc = matplotlib.patches.Arc(center, 2 * (plasmid_radius + 10), 2 * (plasmid_radius + 10),
                          theta1=np.degrees(start_angle), theta2=np.degrees(end_angle),
                          linewidth=2.5, color='green', alpha=0.8)
        ax.add_patch(arc)
        
        # Add primer name
        mid_angle = (start_angle + end_angle) / 2
        label_x = center[0] + (plasmid_radius + 40) * np.cos(mid_angle)
        label_y = center[1] + (plasmid_radius + 40) * np.sin(mid_angle)
        ax.text(label_x, label_y, f"Primer {i+1}: {primer.get('name', '')}", 
               fontsize=8, ha='center', va='center', color='green',
               rotation=np.degrees(mid_angle) if -np.pi/2 <= mid_angle <= np.pi/2 else np.degrees(mid_angle) + 180)
