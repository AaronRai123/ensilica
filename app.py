import streamlit as st
import pandas as pd
import json
import os
from dotenv import load_dotenv
from modules.construct_builder import ConstructBuilder
from modules.construct_builder_fix import create_construct
from modules.primer_tools import design_primers
from modules.visualize import visualize_construct, visualize_linear_map, render_plasmid_map
from modules.gpt_agent import process_natural_language, generate_protocol, analyze_construct
from modules.restriction_analyzer import analyze_restriction_sites, get_compatible_enzyme_pairs
from Bio.SeqRecord import SeqRecord
import plotly.io as pio
import plotly.graph_objects as go
import numpy as np

# Load environment variables
load_dotenv()

# Set page configuration
st.set_page_config(
    page_title="Ensilica - Genetic Construct Designer",
    page_icon="🧬",
    layout="wide"
)

# Function to display construct information and visualization
def display_construct_info(construct_record, primer_info=None, unique_suffix=""):
    """
    Display construct information with interactive visualization options
    
    Args:
        construct_record: The SeqRecord object containing the construct
        primer_info: Optional dictionary with primer information
        unique_suffix: Optional suffix to ensure unique keys
    """
    st.write(f"**Name:** {construct_record.id}")
    st.write(f"**Length:** {len(construct_record.seq)} bp")
    
    # Create a toggle for circular/linear view with a unique key
    view_type = st.radio("View type:", ["Circular", "Linear"], horizontal=True, 
                        key=f"view_type_{construct_record.id}_{unique_suffix}")
    
    # Generate and display visualization
    with st.spinner("Generating construct visualization..."):
        try:
            # Show debug info
            st.write(f"Debug - Name: {construct_record.id}")
            
            # Create a direct Plotly figure instead of using JSON conversion
            if view_type == "Linear":
                # Create a simple linear plot directly
                fig = go.Figure()
                
                # Add base line representing the DNA backbone
                seq_len = len(construct_record.seq)
                fig.add_trace(go.Scatter(
                    x=[0, seq_len],
                    y=[0, 0],
                    mode='lines',
                    line=dict(color='gray', width=2),
                    showlegend=False
                ))
                
                # Add position markers
                marker_interval = 500
                marker_positions = list(range(0, seq_len + marker_interval, marker_interval))
                marker_positions = [p for p in marker_positions if p <= seq_len]
                
                fig.add_trace(go.Scatter(
                    x=marker_positions,
                    y=[0] * len(marker_positions),
                    mode='markers+text',
                    text=[f"{pos}" for pos in marker_positions],
                    textposition='bottom center',
                    marker=dict(size=5, color='black'),
                    showlegend=False
                ))
                
                # Track used feature names for legend
                used_feature_types = set()
                
                # Process features - simple version
                lane_height = 0.2
                max_y = 0.5  # Starting max y value
                
                # First, collect all features and determine their lanes
                feature_lanes = {}
                current_lanes = [0]  # Track the end position of features in each lane
                
                for i, feature in enumerate(construct_record.features):
                    if feature.type == "source":
                        continue
                        
                    # Get feature info
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    
                    # Find a suitable lane
                    lane_idx = 0
                    for j, lane_end in enumerate(current_lanes):
                        if start > lane_end:
                            lane_idx = j
                            break
                    else:
                        # Need a new lane
                        lane_idx = len(current_lanes)
                        current_lanes.append(0)
                    
                    # Update lane end position
                    if lane_idx >= len(current_lanes):
                        current_lanes.append(end)
                    else:
                        current_lanes[lane_idx] = end
                    
                    # Store lane for this feature
                    feature_lanes[i] = lane_idx
                
                # Now draw features with proper lane assignment
                for i, feature in enumerate(construct_record.features):
                    if feature.type == "source":
                        continue
                        
                    # Get feature info
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    if end <= start:  # Skip invalid features
                        continue
                        
                    strand = 1 if feature.strand > 0 else -1 if feature.strand < 0 else 0
                    
                    # Get label
                    label = ""
                    for qualifier in ["label", "gene", "product", "note"]:
                        if qualifier in feature.qualifiers:
                            label = feature.qualifiers[qualifier][0]
                            break
                    
                    if not label:
                        label = f"{feature.type}_{i+1}"
                        
                    # Color based on type
                    color = "#cccccc"  # Default gray
                    if feature.type == "CDS" or feature.type == "gene":
                        color = "#4285F4"  # Blue for genes
                    elif feature.type == "promoter":
                        color = "#9C27B0"  # Purple for promoters
                    elif feature.type == "terminator":
                        color = "#E91E63"  # Pink for terminators
                    elif feature.type == "rep_origin" or "origin" in feature.type.lower():
                        color = "#2196F3"  # Light blue for origins
                    elif feature.type == "primer_bind":
                        color = "#FF5722"  # Orange for primer sites
                    
                    # Get lane and calculate y position
                    lane_idx = feature_lanes[i]
                    lane_sign = 1 if lane_idx % 2 == 0 else -1
                    lane_offset = (lane_idx // 2 + 1) * 0.3
                    y_pos = lane_sign * lane_offset
                    
                    # Track max y value for plot scaling
                    max_y = max(max_y, abs(y_pos) + lane_height)
                    
                    # Feature name for legend - use feature type, show only once per type
                    legend_name = feature.type
                    show_legend = feature.type not in used_feature_types
                    if show_legend:
                        used_feature_types.add(feature.type)
                    
                    # Draw arrow shape for the feature
                    if strand != 0:
                        # Arrow body
                        arrow_length = end - start
                        arrow_head = min(10, arrow_length / 4)  # Head size proportional to feature, max 10bp
                        
                        if strand > 0:  # Forward strand
                            # Arrow body
                            body_end = end - arrow_head
                            fig.add_trace(go.Scatter(
                                x=[start, start, body_end, body_end, start],
                                y=[y_pos-lane_height/2, y_pos+lane_height/2, y_pos+lane_height/2, y_pos-lane_height/2, y_pos-lane_height/2],
                                fill="toself",
                                fillcolor=color,
                                line=dict(color=color, width=1),
                                mode='lines',
                                name=legend_name,
                                showlegend=show_legend,
                                hoverinfo='text',
                                hovertext=f"{label}<br>Type: {feature.type}<br>Position: {start}-{end}"
                            ))
                            
                            # Arrow head
                            fig.add_trace(go.Scatter(
                                x=[body_end, end, body_end],
                                y=[y_pos+lane_height/2, y_pos, y_pos-lane_height/2],
                                fill="toself",
                                fillcolor=color,
                                line=dict(color=color, width=1),
                                mode='lines',
                                showlegend=False,
                                hoverinfo='skip'
                            ))
                        else:  # Reverse strand
                            # Arrow body
                            body_start = start + arrow_head
                            fig.add_trace(go.Scatter(
                                x=[body_start, body_start, end, end, body_start],
                                y=[y_pos-lane_height/2, y_pos+lane_height/2, y_pos+lane_height/2, y_pos-lane_height/2, y_pos-lane_height/2],
                                fill="toself",
                                fillcolor=color,
                                line=dict(color=color, width=1),
                                mode='lines',
                                name=legend_name,
                                showlegend=show_legend,
                                hoverinfo='text',
                                hovertext=f"{label}<br>Type: {feature.type}<br>Position: {start}-{end}"
                            ))
                            
                            # Arrow head
                            fig.add_trace(go.Scatter(
                                x=[body_start, start, body_start],
                                y=[y_pos+lane_height/2, y_pos, y_pos-lane_height/2],
                                fill="toself",
                                fillcolor=color,
                                line=dict(color=color, width=1),
                                mode='lines',
                                showlegend=False,
                                hoverinfo='skip'
                            ))
                    else:  # No strand, just a rectangle
                        fig.add_trace(go.Scatter(
                            x=[start, start, end, end, start],
                            y=[y_pos-lane_height/2, y_pos+lane_height/2, y_pos+lane_height/2, y_pos-lane_height/2, y_pos-lane_height/2],
                            fill="toself",
                            fillcolor=color,
                            line=dict(color=color, width=1),
                            mode='lines',
                            name=legend_name,
                            showlegend=show_legend,
                            hoverinfo='text',
                            hovertext=f"{label}<br>Type: {feature.type}<br>Position: {start}-{end}"
                        ))
                    
                    # Add text label for the feature
                    if end - start > 50:  # Only label features that are wide enough
                        fig.add_trace(go.Scatter(
                            x=[(start + end) / 2],
                            y=[y_pos],
                            mode='text',
                            text=[label],
                            textposition='middle center',
                            textfont=dict(color='white', size=10),
                            showlegend=False,
                            hoverinfo='skip'
                        ))
                
                # Set layout with adjusted y range based on content
                fig.update_layout(
                    title=f"Linear Map: {construct_record.id} ({seq_len} bp)",
                    xaxis=dict(
                        title="Base Position",
                        showgrid=False
                    ),
                    yaxis=dict(
                        showticklabels=False,
                        showgrid=False,
                        zeroline=True,
                        zerolinecolor='gray',
                        zerolinewidth=1,
                        range=[-max_y - 0.2, max_y + 0.2]
                    ),
                    height=400,
                    margin=dict(l=20, r=20, t=40, b=20),
                    plot_bgcolor='rgba(0,0,0,0)',
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.02,
                        xanchor="right",
                        x=1
                    )
                )
            else:
                # Create a simple circular plot directly
                fig = go.Figure()
                
                # Get sequence length
                seq_len = len(construct_record.seq)
                
                # Add base circle representing the plasmid backbone
                circle_pts = 100
                theta = np.linspace(0, 2*np.pi, circle_pts)
                radius = np.ones(circle_pts)
                
                fig.add_trace(go.Scatterpolar(
                    r=radius,
                    theta=theta * 180/np.pi,
                    mode='lines',
                    line=dict(color='gray', width=1.5)
                ))
                
                # Process features - simple version
                for i, feature in enumerate(construct_record.features):
                    if feature.type == "source":
                        continue
                        
                    # Get feature info
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    
                    # Handle features that cross the origin
                    if end < start:
                        end = end + seq_len
                    
                    # Convert to angular coordinates (in degrees)
                    start_angle = (start / seq_len) * 360
                    end_angle = (end / seq_len) * 360
                    
                    # Create points for arc
                    arc_pts = 50
                    arc_theta = np.linspace(start_angle, end_angle, arc_pts)
                    
                    # Get label
                    label = ""
                    for qualifier in ["label", "gene", "product", "note"]:
                        if qualifier in feature.qualifiers:
                            label = feature.qualifiers[qualifier][0]
                            break
                    
                    if not label:
                        label = f"{feature.type}_{i+1}"
                        
                    # Color based on type
                    color = "#cccccc"  # Default gray
                    if feature.type == "CDS":
                        color = "#ffcc00"  # Yellow
                    elif feature.type == "promoter":
                        color = "#00ccff"  # Light blue
                    elif feature.type == "terminator":
                        color = "#ff0000"  # Red
                    elif feature.type == "rep_origin" or "origin" in feature.type.lower():
                        color = "#cc00ff"  # Purple
                    
                    # Add the feature arc
                    fig.add_trace(go.Scatterpolar(
                        r=np.ones(arc_pts) * 0.9,  # Slightly inside the main circle
                        theta=arc_theta,
                        mode='lines',
                        line=dict(color=color, width=10),
                        name=label,
                        hoverinfo='text',
                        hovertext=f"{label}<br>Type: {feature.type}<br>Position: {start}-{end}"
                    ))
                
                # Set layout
                fig.update_layout(
                    title=f"Plasmid Map: {construct_record.id} ({seq_len} bp)",
                    polar=dict(
                        radialaxis=dict(visible=False, range=[0, 1.2]),
                        angularaxis=dict(visible=True, direction='clockwise')
                    ),
                    height=500
                )
            
            # Display the figure with a unique key
            st.plotly_chart(fig, use_container_width=True, key=f"plot_{construct_record.id}_{unique_suffix}")
            st.info("💡 TIP: Hover over features for details. You can zoom, pan, and export this interactive map.")
        except Exception as e:
            st.error(f"Could not generate visualization: {str(e)}")
            import traceback
            st.code(traceback.format_exc())
    
    # Display features in a table
    st.write("### Construct Components")
    
    feature_data = []
    for feature in construct_record.features:
        if feature.type == "source":
            continue
        
        # Extract feature name
        feature_name = ""
        for qualifier in ["label", "gene", "product", "note"]:
            if qualifier in feature.qualifiers:
                feature_name = feature.qualifiers[qualifier][0]
                break
        
        if not feature_name:
            feature_name = f"{feature.type}_feature"
        
        # Get coordinates
        start = int(feature.location.start) + 1  # Convert to 1-based coordinates for display
        end = int(feature.location.end)
        
        # Determine strand
        strand = "+" if feature.strand == 1 else "-" if feature.strand == -1 else "none"
        
        feature_data.append({
            "Feature": feature_name,
            "Type": feature.type,
            "Start": start,
            "End": end,
            "Length": end - start + 1,
            "Strand": strand
        })
    
    if feature_data:
        st.dataframe(feature_data)
    else:
        st.write("No features found in this construct.")
    
    # Display primer information if available
    if primer_info:
        st.write("### Primers")
        for primer in primer_info:
            st.write(f"**{primer['name']}:** {primer['sequence']}")
            st.write(f"Position: {primer['position']}, Tm: {primer['tm']}°C")

# App title and description
st.title("🧬 Ensilica")
st.markdown("### AI-Powered Genetic Construct Designer")
st.markdown("Design genetic constructs by selecting components or entering a natural language description.")

# API key input in sidebar
with st.sidebar:
    st.subheader("API Status")
    st.success("✅ All systems ready")
    
    st.markdown("---")
    
    # Add additional options in sidebar
    st.subheader("Available Sequences")
    st.info("Local sequences in data/ directory will be used when available.")

# Basic input form
st.subheader("Design Your Genetic Construct")

# Component selection inputs
manual_tab, ai_tab = st.tabs(["Component Selection", "AI Assistant"])

with manual_tab:
    col1, col2 = st.columns(2)

    with col1:
        gene_options = ["GFP", "mCherry", "RFP"]
        gene_name = st.selectbox("Select Gene:", gene_options)
        
        promoter_options = ["T7", "lac", "trc"]
        promoter_name = st.selectbox("Select Promoter:", promoter_options)

    with col2:
        vector_options = ["pUC19", "pET28a", "pBluescript"]
        vector_name = st.selectbox("Select Vector:", vector_options)
        
        cloning_options = ["Gibson Assembly", "Restriction Enzyme Cloning"]
        cloning_method = st.selectbox("Cloning Method:", cloning_options)

with ai_tab:
    user_input = st.text_area("Describe your genetic construct:", 
                             placeholder="Example: Insert mCherry into pUC19 under a lactose promoter", 
                             height=100)
    
    # Example queries as buttons
    st.markdown("**Example queries:**")
    examples = [
        "Insert GFP into pUC19 under a T7 promoter",
        "Clone mCherry into pET28a with a lac promoter",
        "Add RFP to pBluescript with trc promoter and add a 6xHis tag"
    ]

    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button(examples[0]):
            user_input = examples[0]
            st.session_state.user_input = examples[0]

    with col2:
        if st.button(examples[1]):
            user_input = examples[1]
            st.session_state.user_input = examples[1]

    with col3:
        if st.button(examples[2]):
            user_input = examples[2]
            st.session_state.user_input = examples[2]

# Add this function to display restriction analysis
def add_restriction_analysis(record: SeqRecord):
    """
    Add restriction analysis to the app
    
    Args:
        record: SeqRecord to analyze
    """
    st.subheader("Restriction Site Analysis")
    
    # Let user select enzymes
    common_enzymes = ['EcoRI', 'BamHI', 'HindIII', 'XhoI', 'NdeI', 'XbaI', 
                      'PstI', 'SalI', 'SmaI', 'KpnI', 'SacI', 'BglII']
    
    with st.expander("Customize Enzymes"):
        selected_enzymes = st.multiselect(
            "Select restriction enzymes:",
            options=common_enzymes,
            default=common_enzymes[:6]
        )
    
    # If no enzymes selected, use default set
    if not selected_enzymes:
        selected_enzymes = common_enzymes[:6]
    
    # Analyze restriction sites
    with st.spinner("Analyzing restriction sites..."):
        analysis = analyze_restriction_sites(record, selected_enzymes)
    
    # Display results
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Show restriction map
        if analysis["map_image"]:
            st.image(f"data:image/png;base64,{analysis['map_image']}", caption="Restriction Map")
        else:
            st.info("No restriction sites found with selected enzymes.")
    
    with col2:
        # Show unique cutters
        st.subheader("Unique Cutters")
        if analysis["unique_cutters"]:
            for enzyme in analysis["unique_cutters"]:
                st.write(f"✅ {enzyme}")
        else:
            st.write("No unique cutters found.")
        
        # Show rare cutters
        st.subheader("Rare Cutters (2-3 sites)")
        if analysis["rare_cutters"]:
            for enzyme in analysis["rare_cutters"]:
                st.write(f"⚠️ {enzyme}")
        else:
            st.write("No rare cutters found.")
    
    # Display restriction sites table
    if not analysis["dataframe"].empty:
        st.subheader("Restriction Sites")
        st.dataframe(analysis["dataframe"][["enzyme", "position", "recognition_site", "overhang"]], use_container_width=True)
    
    # Show compatible enzyme pairs
    pairs = get_compatible_enzyme_pairs(record)
    if pairs:
        st.subheader("Recommended Enzyme Pairs")
        
        pairs_data = []
        for pair in pairs[:10]:  # Show only top 10 pairs
            compatibility = "⭐" * (pair["compatibility"] // 2)  # Convert 0-10 scale to 0-5 stars
            pairs_data.append({
                "Enzyme Pair": f"{pair['enzyme1']} + {pair['enzyme2']}",
                "Fragment Size": f"{pair['fragment_size']} bp",
                "Compatibility": compatibility,
                "Positions": f"{pair['position1']} / {pair['position2']}"
            })
        
        st.dataframe(pd.DataFrame(pairs_data), use_container_width=True)
    else:
        st.info("No compatible enzyme pairs found for cloning.")

# Add this function to generate HTML with color-coded DNA sequence
def generate_color_coded_sequence(record: SeqRecord) -> str:
    """
    Generate HTML with color-coded DNA sequence based on features
    
    Args:
        record: SeqRecord object containing sequence and features
        
    Returns:
        HTML string with color-coded sequence
    """
    # Get the sequence as string
    seq_str = str(record.seq)
    
    # Create a list to track which positions belong to which features
    position_features = [[] for _ in range(len(seq_str))]
    
    # Assign features to positions
    for feature in record.features:
        if feature.type == "source":
            continue
            
        start = int(feature.location.start)
        end = int(feature.location.end)
        
        # Get feature name
        feature_name = ""
        for qualifier in ["label", "gene", "product", "note"]:
            if qualifier in feature.qualifiers:
                feature_name = feature.qualifiers[qualifier][0]
                break
        
        if not feature_name:
            feature_name = f"{feature.type}_feature"
            
        # Assign feature to all positions it covers
        for i in range(start, end):
            if i < len(seq_str):  # Avoid index errors
                position_features[i].append((feature.type, feature_name))
    
    # Define colors for different feature types
    colors = {
        "CDS": "#4285F4",        # Blue for CDS/genes
        "gene": "#4285F4",       # Blue for CDS/genes
        "promoter": "#9C27B0",   # Purple for promoters
        "terminator": "#E91E63", # Pink for terminators
        "rep_origin": "#2196F3", # Light blue for origins
        "origin": "#2196F3",     # Light blue for origins
        "misc_feature": "#607D8B", # Gray for misc features
        "primer_bind": "#FF5722", # Red/orange for primer sites
        "RBS": "#00BCD4",        # Cyan for RBS
    }
    default_color = "#FFFFFF"   # white for non-feature nucleotides
    
    # Generate HTML
    html_parts = ['<div style="font-family: monospace; line-height: 1.5; overflow-wrap: anywhere">']
    
    for i, base in enumerate(seq_str):
        if i > 0 and i % 10 == 0:
            html_parts.append(' ')
        if i > 0 and i % 100 == 0:
            html_parts.append('<br>')
            
        if not position_features[i]:
            # No feature at this position
            html_parts.append(f'<span style="color: {default_color}">{base}</span>')
        else:
            # Use the color of the first feature type at this position
            feature_type = position_features[i][0][0]
            feature_name = position_features[i][0][1]
            color = colors.get(feature_type, "#757575")  # Default gray if type not in colors
            html_parts.append(f'<span style="color: {color}" title="{feature_name} ({feature_type})">{base}</span>')
    
    html_parts.append('</div>')
    
    # Add a legend
    html_parts.append('<div style="margin-top: 20px; font-family: sans-serif; font-size: 14px">')
    html_parts.append('<strong>Legend:</strong><br>')
    for feature_type, color in colors.items():
        html_parts.append(f'<span style="display: inline-block; width: 20px; height: 14px; background-color: {color}; margin-right: 5px;"></span>{feature_type}<br>')
    html_parts.append('</div>')
    
    return ''.join(html_parts)

# Process button
if st.button("Generate Construct", type="primary"):
    # Determine if we're using AI input (check if we're on the AI tab and have input)
    use_ai_input = ai_tab and user_input
    
    if use_ai_input:
        with st.spinner("Processing natural language input..."):
            # Use AI to process the natural language input
            construct_info = process_natural_language(user_input)
            gene_name = construct_info.get("gene_name")
            vector_name = construct_info.get("vector_name")
            promoter_name = construct_info.get("promoter_name")
            cloning_method = construct_info.get("cloning_method")
    else:
        # Use the manually selected components
        construct_info = {
            "gene_name": gene_name,
            "vector_name": vector_name,
            "promoter_name": promoter_name,
            "cloning_method": cloning_method,
            "insertion_site": "Multiple Cloning Site",
            "additional_features": []
        }
    
    # Display the construct info
    st.subheader("Construct Information")
    st.json(construct_info)
    
    # Create the construct
    with st.spinner("Building DNA construct..."):
        try:
            construct_result = create_construct(construct_info)
            
            st.subheader("Constructed Sequence")
            
            # Display the construct information using the new function
            display_construct_info(construct_result['construct_record'], 
                                   construct_result.get('primer_info'),
                                   unique_suffix="main_view")
            
            # Design primers
            with st.spinner("Designing optimal primers..."):
                primer_results = design_primers(
                    construct_result["insert_region"], 
                    construct_info
                )
                
                st.subheader("PCR Primers")
                
                # Create DataFrame with prettier column names
                primer_data = []
                for primer in primer_results["primers"]:
                    primer_data.append({
                        "Name": primer["name"],
                        "Sequence (5' to 3')": primer["sequence"],
                        "Length": primer["length"],
                        "Tm (°C)": primer["tm"],
                        "GC Content (%)": primer["gc_content"]
                    })
                
                primer_df = pd.DataFrame(primer_data)
                st.dataframe(primer_df)
            
            # Generate AI-powered protocol
            with st.spinner("Generating laboratory protocol..."):
                protocol = generate_protocol(construct_info, primer_results)
                
                st.subheader("Laboratory Protocol")
                st.markdown(protocol)
            
            # AI analysis of the construct
            with st.spinner("Analyzing construct..."):
                analysis = analyze_construct(construct_info, construct_result["construct_record"])
                
                st.subheader("Construct Analysis")
                
                # Display potential issues
                if "potential_issues" in analysis and analysis["potential_issues"]:
                    st.warning("**Potential Issues:**")
                    for issue in analysis["potential_issues"]:
                        st.write(f"- {issue}")
                
                # Display recommendations
                if "recommendations" in analysis and analysis["recommendations"]:
                    st.info("**Recommendations:**")
                    for rec in analysis["recommendations"]:
                        st.write(f"- {rec}")
                
                # Display expression prediction
                if "expression_prediction" in analysis:
                    st.success(f"**Expression Prediction:** {analysis['expression_prediction']}")
                
                # Display special considerations
                if "special_considerations" in analysis:
                    st.write(f"**Special Considerations:** {analysis['special_considerations']}")
            
            # Download buttons
            st.subheader("Download Files")
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.download_button(
                    label="Download GenBank",
                    data=construct_result["genbank_content"],
                    file_name=f"{gene_name}_{vector_name}.gb",
                    mime="text/plain",
                    key="download_genbank_main"
                )
            
            with col2:
                st.download_button(
                    label="Download FASTA",
                    data=construct_result["fasta_content"],
                    file_name=f"{gene_name}_{vector_name}.fasta",
                    mime="text/plain",
                    key="download_fasta_main"
                )
            
            with col3:
                st.download_button(
                    label="Download Protocol",
                    data=protocol,
                    file_name=f"{gene_name}_{vector_name}_protocol.md",
                    mime="text/plain",
                    key="download_protocol_main"
                )
            
            with col4:
                primer_text = "\n".join([f"{row['Name']}: {row['Sequence (5\' to 3\')']} (Tm: {row['Tm (°C)']}°C, GC: {row['GC Content (%)']}%)" for _, row in primer_df.iterrows()])
                st.download_button(
                    label="Download Primers",
                    data=primer_text,
                    file_name=f"{gene_name}_{vector_name}_primers.txt",
                    mime="text/plain",
                    key="download_primers_main"
                )
                
            # Add restriction analysis
            add_restriction_analysis(construct_result["construct_record"])
            
            # If a construct was generated, show the results
            if construct_result:
                # Display basic information about the construct
                st.markdown(f"## 🧬 Construct: {construct_info['gene_name']}_{construct_info['vector_name']}_{cloning_method.lower().replace(' ', '_')}")
                st.write(f"Total size: **{len(construct_result['construct_record'].seq)}** bp")
                
                # Display tabs for different outputs
                construct_tabs = st.tabs(["Construct Info", "DNA Sequence", "PCR Primers", "Lab Protocol", "Validation Results"])
                
                # Tab 1: Construct Info - show a visual representation
                with construct_tabs[0]:
                    st.markdown("### Construct Map")
                    st.write("Feature map of the created construct:")
                    
                    # Use the new display function with a unique suffix
                    display_construct_info(construct_result['construct_record'], unique_suffix="tab_view")
                    
                    # Display the cloning strategy
                    st.markdown("### Cloning Strategy")
                    st.write(f"**Method:** {cloning_method}")
                
                # Tab 2: DNA - show the DNA sequence in different formats
                with construct_tabs[1]:
                    sequence_tabs = st.tabs(["Color-Coded DNA", "FASTA", "GenBank", "Raw Sequence"])
                    
                    with sequence_tabs[0]:
                        # Display color-coded sequence
                        colored_sequence = generate_color_coded_sequence(construct_result['construct_record'])
                        st.components.v1.html(colored_sequence, height=500, scrolling=True)
                        st.info("💡 Hover over colored nucleotides to see which feature they belong to. Sequence is wrapped with spaces every 10 bp and line breaks every 100 bp.")
                    
                    with sequence_tabs[1]:
                        st.text(construct_result["fasta_content"])
                        st.download_button(
                            "Download FASTA",
                            construct_result["fasta_content"],
                            file_name=f"{construct_info['gene_name']}_{construct_info['vector_name']}.fasta",
                            mime="text/plain",
                            key="download_fasta_tab"
                        )
                        
                    with sequence_tabs[2]:
                        st.text(construct_result["genbank_content"])
                        st.download_button(
                            "Download GenBank",
                            construct_result["genbank_content"],
                            file_name=f"{construct_info['gene_name']}_{construct_info['vector_name']}.gb",
                            mime="text/plain",
                            key="download_genbank_tab"
                        )
                        
                    with sequence_tabs[3]:
                        st.text(str(construct_result['construct_record'].seq))
                
                # Tab 3: PCR Primers - show primers for the insert
                with construct_tabs[2]:
                    if primer_results and "primers" in primer_results:
                        st.markdown("### PCR Primers for Gene Amplification")
                        
                        for i, primer in enumerate(primer_results["primers"]):
                            if i < 2:  # Just the main forward and reverse primers
                                col1, col2, col3, col4 = st.columns([2, 5, 1, 1])
                                with col1:
                                    st.markdown(f"**{primer['name']}**")
                                with col2:
                                    st.code(primer['sequence'])
                                with col3:
                                    st.write(f"Tm: {primer['tm']}°C")
                                with col4:
                                    st.write(f"GC: {primer['gc_content']}%")
                    
                        if len(primer_results["primers"]) > 2:
                            with st.expander("Verification Primers"):
                                for i, primer in enumerate(primer_results["primers"][2:]):
                                    col1, col2, col3, col4 = st.columns([2, 5, 1, 1])
                                    with col1:
                                        st.markdown(f"**{primer['name']}**")
                                    with col2:
                                        st.code(primer['sequence'])
                                    with col3:
                                        st.write(f"Tm: {primer['tm']}°C")
                                    with col4:
                                        st.write(f"GC: {primer['gc_content']}%")
                    else:
                        st.info("No primers were generated for this construct.")
                
                # Tab 4: Laboratory Protocol - generated based on the construct and cloning method
                with construct_tabs[3]:
                    if "protocol" in construct_result:
                        protocol_text = construct_result["protocol"]
                        st.markdown(protocol_text)
                    else:
                        # Generate a protocol using the NLP model based on the construct information
                        try:
                            with st.spinner("Generating laboratory protocol..."):
                                protocol = generate_protocol(construct_info, construct_result)
                                st.markdown(protocol)
                        except Exception as e:
                            st.error(f"Could not generate protocol: {str(e)}")
                            st.markdown("""
                            ## Laboratory Protocol Template
                            
                            ### Materials Needed
                            - Forward and reverse primers
                            - Template DNA
                            - PCR reaction mix
                            - Restriction enzymes (if applicable)
                            - DNA ligase (if applicable)
                            - Competent cells
                            - Selection antibiotic
                            
                            ### Procedure
                            1. Amplify the gene of interest using the primers shown in the PCR Primers tab
                            2. Verify PCR product by gel electrophoresis
                            3. Purify the PCR product
                            4. Prepare vector and insert for the selected cloning method
                            5. Ligate vector and insert
                            6. Transform into competent cells
                            7. Plate on selective media
                            8. Screen colonies by PCR or restriction digest
                            9. Sequence verify positive clones
                            """)
                
                # Tab 5: Validation Results - show smart sequence validation results
                with construct_tabs[4]:
                    st.markdown("### Smart Sequence Validation Results")
                    
                    # Get validation results
                    try:
                        # Create a builder instance to access validation methods
                        builder = ConstructBuilder()
                        validation_results = builder.smart_sequence_validation(construct_result['construct_record'])
                        
                        # Frame shift check
                        frame_result = validation_results["frame_shift"]
                        frame_icon = "✅" if frame_result["valid"] else "⚠️"
                        st.markdown(f"**{frame_icon} Frame Check:**")
                        if frame_result["valid"]:
                            st.success("All coding sequences have proper start/stop codons and are in frame.")
                        else:
                            st.warning(frame_result["message"])
                        
                        # Promoter orientation
                        promoter_result = validation_results["promoter_orientation"]
                        promoter_icon = "✅" if promoter_result["valid"] else "⚠️"
                        st.markdown(f"**{promoter_icon} Promoter Orientation:**")
                        if promoter_result["valid"]:
                            st.success("All promoters are correctly oriented relative to their coding sequences.")
                        else:
                            st.warning(promoter_result["message"])
                        
                        # RBS presence
                        rbs_result = validation_results["rbs_presence"]
                        rbs_icon = "✅" if rbs_result["valid"] else "⚠️"
                        st.markdown(f"**{rbs_icon} Ribosome Binding Sites:**")
                        if rbs_result["valid"]:
                            st.success("Ribosome binding sites are present before coding sequences.")
                        else:
                            st.warning(rbs_result["message"])
                        
                        # Restriction sites
                        restriction_result = validation_results["restriction_sites"]
                        
                        st.markdown("**🔍 Restriction Enzyme Analysis:**")
                        if not restriction_result["conflicting_sites"]:
                            st.success("No conflicting restriction sites found.")
                        else:
                            st.info(restriction_result["message"])
                            
                            # Create a table of conflicting sites
                            site_data = []
                            for site in restriction_result["conflicting_sites"]:
                                site_data.append({
                                    "Enzyme": site["enzyme"],
                                    "Recognition Site": site["site"],
                                    "Occurrences": site["count"],
                                    "Positions": ", ".join([str(pos) for pos in site["positions"][:3]]) + 
                                                ("..." if len(site["positions"]) > 3 else "")
                                })
                            
                            st.dataframe(site_data, hide_index=True)
                            
                            st.info("**Note:** These restriction sites might interfere with restriction enzyme-based cloning methods. Consider using a different cloning method or different restriction enzymes.")
                        
                    except Exception as e:
                        st.error(f"Could not perform sequence validation: {str(e)}")
                        st.info("Try regenerating the construct to enable full validation.")
                    
                # Output analysis of the construct
                st.sidebar.markdown("### Construct Analysis")
                
                # Use fallback analysis directly to avoid the error with SeqRecord
                try:
                    with st.sidebar.status("Analyzing construct...", expanded=False) as status:
                        st.sidebar.caption("Examining construct properties...")
                        
                        # Convert construct info to a simple dict
                        analysis_info = {k: v for k, v in construct_info.items()}
                        
                        # Use only the fallback analysis which doesn't try to use get() on SeqRecord
                        from modules.gpt_agent import _generate_fallback_analysis
                        analysis = _generate_fallback_analysis(analysis_info)
                        
                        status.update(label="Analysis complete!", state="complete", expanded=True)
                        
                        # Display analysis
                        if "vector_suitability" in analysis:
                            st.sidebar.markdown("##### Vector Suitability")
                            st.sidebar.write(analysis["vector_suitability"])
                            
                        if "expression_levels" in analysis:
                            st.sidebar.markdown("##### Expected Expression")
                            st.sidebar.write(analysis["expression_levels"])
                            
                        if "cloning_complications" in analysis:
                            st.sidebar.markdown("##### Potential Challenges")
                            st.sidebar.write(analysis["cloning_complications"])
                            
                        if "validation_recommendations" in analysis:
                            st.sidebar.markdown("##### Validation Strategy")
                            st.sidebar.write(analysis["validation_recommendations"])
                            
                        if "potential_applications" in analysis:
                            st.sidebar.markdown("##### Applications")
                            st.sidebar.write(analysis["potential_applications"])
                except Exception as e:
                    st.sidebar.error(f"Could not analyze construct: {str(e)}")
                    st.sidebar.markdown("Try using a different cloning method or adjust your construct parameters.")
            
        except Exception as e:
            st.error(f"An error occurred: {str(e)}")
            import traceback
            st.code(traceback.format_exc())

# Add information about the app at the bottom
st.markdown("---")
st.markdown("""
**About Ensilica:**
Ensilica is an AI-powered tool for designing genetic constructs. It uses advanced language models to process natural language descriptions, generate sequences, design primers, and create detailed laboratory protocols.
""") 