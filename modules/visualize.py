"""
Visualization module for DNA constructs - simplified for compatibility
"""
import base64
import io
import math
from typing import Dict, Any, Optional
from Bio.SeqRecord import SeqRecord

def visualize_construct(construct_record: SeqRecord, insert_region: Optional[Dict[str, Any]] = None) -> str:
    """
    Generate a simplified text-based visualization of the DNA construct.
    
    Args:
        construct_record: SeqRecord object containing the construct
        insert_region: Optional dictionary with information about the inserted region
        
    Returns:
        Base64-encoded SVG image
    """
    try:
        # Generate an SVG visualization
        return create_svg_visualization(construct_record, insert_region)
    except Exception as e:
        print(f"Error generating visualization: {str(e)}")
        # Return a placeholder SVG
        return create_placeholder_svg(construct_record)

def create_svg_visualization(construct_record: SeqRecord, insert_region: Optional[Dict[str, Any]] = None) -> str:
    """Create an SVG visualization of the plasmid."""
    sequence_length = len(construct_record.seq)
    
    # SVG settings
    width, height = 500, 500
    center_x, center_y = width / 2, height / 2
    radius = 180
    
    # Start SVG
    svg = f"""
    <svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
        <circle cx="{center_x}" cy="{center_y}" r="{radius}" fill="white" stroke="black" stroke-width="2"/>
        <text x="{center_x}" y="{center_y - radius - 20}" text-anchor="middle" font-size="16" font-weight="bold">{construct_record.id} ({sequence_length} bp)</text>
    """
    
    # Add features
    feature_colors = {
        'CDS': 'lightblue',
        'promoter': 'lightgreen',
        'terminator': 'salmon',
        'rep_origin': 'gold',
        'misc_feature': 'lightgray'
    }
    
    for i, feature in enumerate(construct_record.features):
        feature_type = feature.type
        feature_color = feature_colors.get(feature_type, 'lightgray')
        
        # Check if this is the inserted gene
        is_insert = False
        if insert_region and feature.location.start <= insert_region['end'] and feature.location.end >= insert_region['start']:
            feature_color = 'tomato'  # Highlight inserted region
            is_insert = True
        
        # Get feature name
        feature_name = ""
        if "label" in feature.qualifiers:
            feature_name = feature.qualifiers["label"][0]
        elif "product" in feature.qualifiers:
            feature_name = feature.qualifiers["product"][0]
        else:
            feature_name = f"{feature_type}_{i+1}"
        
        # Calculate angles for the feature
        start_bp = int(feature.location.start)
        end_bp = int(feature.location.end)
        
        start_angle = 2 * math.pi * start_bp / sequence_length
        end_angle = 2 * math.pi * end_bp / sequence_length
        
        # Draw arc/wedge for the feature
        feature_width = 30
        inner_radius = radius - feature_width
        
        # Calculate path for the arc
        start_x = center_x + radius * math.cos(start_angle)
        start_y = center_y + radius * math.sin(start_angle)
        end_x = center_x + radius * math.cos(end_angle)
        end_y = center_y + radius * math.sin(end_angle)
        inner_start_x = center_x + inner_radius * math.cos(start_angle)
        inner_start_y = center_y + inner_radius * math.sin(start_angle)
        inner_end_x = center_x + inner_radius * math.cos(end_angle)
        inner_end_y = center_y + inner_radius * math.sin(end_angle)
        
        # Determine if arc is more than half the circle
        large_arc_flag = 1 if end_angle - start_angle > math.pi else 0
        
        # Create path
        path = f"""
        <path d="M {start_x} {start_y} 
                A {radius} {radius} 0 {large_arc_flag} 1 {end_x} {end_y}
                L {inner_end_x} {inner_end_y}
                A {inner_radius} {inner_radius} 0 {large_arc_flag} 0 {inner_start_x} {inner_start_y}
                Z" 
            fill="{feature_color}" stroke="black" stroke-width="1" opacity="0.8"/>
        """
        svg += path
        
        # Add label for the feature
        label_angle = (start_angle + end_angle) / 2
        label_radius = radius + 20
        label_x = center_x + label_radius * math.cos(label_angle)
        label_y = center_y + label_radius * math.sin(label_angle)
        
        # Adjust text rotation for readability
        rotation = label_angle * 180 / math.pi
        if rotation > 90 and rotation < 270:
            rotation += 180
        
        svg += f"""
        <text x="{label_x}" y="{label_y}" 
            text-anchor="middle" 
            transform="rotate({rotation} {label_x} {label_y})"
            font-size="12">{feature_name}</text>
        """
    
    # Close SVG
    svg += "</svg>"
    
    # Convert to base64
    svg_bytes = svg.encode('utf-8')
    return base64.b64encode(svg_bytes).decode('utf-8')

def create_placeholder_svg(construct_record: SeqRecord) -> str:
    """Create a simple placeholder SVG when visualization fails."""
    # Simple SVG
    svg = f"""
    <svg width="400" height="300" xmlns="http://www.w3.org/2000/svg">
        <rect width="400" height="300" fill="white" stroke="black" stroke-width="2"/>
        <text x="200" y="150" text-anchor="middle" font-size="16">{construct_record.id}</text>
        <text x="200" y="180" text-anchor="middle" font-size="14">{len(construct_record.seq)} bp</text>
    </svg>
    """
    
    # Convert to base64
    svg_bytes = svg.encode('utf-8')
    return base64.b64encode(svg_bytes).decode('utf-8') 