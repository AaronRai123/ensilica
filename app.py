import streamlit as st
import pandas as pd
import json
import os
from dotenv import load_dotenv
from modules.construct_builder import create_construct
from modules.primer_tools import design_primers
from modules.visualize import visualize_construct
from modules.gpt_agent import process_natural_language, generate_protocol, analyze_construct
from modules.restriction_analyzer import analyze_restriction_sites, get_compatible_enzyme_pairs
from Bio.SeqRecord import SeqRecord

# Load environment variables
load_dotenv()

# Set page configuration
st.set_page_config(
    page_title="Ensilica - Genetic Construct Designer",
    page_icon="🧬",
    layout="wide"
)

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
            
            # Split into two columns for sequence and visualization
            seq_col, viz_col = st.columns(2)
            
            with seq_col:
                st.text_area("FASTA Sequence (Preview)", construct_result["fasta_content"], height=150)
            
            with viz_col:
                # Generate and display visualization
                viz_image = visualize_construct(
                    construct_result["construct_record"], 
                    construct_result["insert_region"]
                )
                st.image(f"data:image/svg+xml;base64,{viz_image}", caption="Plasmid Map")
            
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
                analysis = analyze_construct(construct_result["construct_record"], construct_info)
                
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
                    mime="text/plain"
                )
            
            with col2:
                st.download_button(
                    label="Download FASTA",
                    data=construct_result["fasta_content"],
                    file_name=f"{gene_name}_{vector_name}.fasta",
                    mime="text/plain"
                )
            
            with col3:
                st.download_button(
                    label="Download Protocol",
                    data=protocol,
                    file_name=f"{gene_name}_{vector_name}_protocol.md",
                    mime="text/plain"
                )
            
            with col4:
                primer_text = "\n".join([f"{row['Name']}: {row['Sequence (5\' to 3\')']} (Tm: {row['Tm (°C)']}°C, GC: {row['GC Content (%)']}%)" for _, row in primer_df.iterrows()])
                st.download_button(
                    label="Download Primers",
                    data=primer_text,
                    file_name=f"{gene_name}_{vector_name}_primers.txt",
                    mime="text/plain"
                )
                
            # Add restriction analysis
            add_restriction_analysis(construct_result["construct_record"])
            
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