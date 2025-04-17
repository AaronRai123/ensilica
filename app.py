import streamlit as st
import pandas as pd
import json
import os
from dotenv import load_dotenv
from modules.construct_builder import create_construct, ConstructBuilder
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
                with st.spinner("Generating construct visualization..."):
                    try:
                        # Pass the insert region to highlight it
                        image_data = visualize_construct(construct_result['construct_record'], construct_result.get('insert_region'))
                        
                        # Use HTML component instead of st.image for better SVG rendering
                        st.markdown(f"""
                        <div style="text-align: center; margin: 10px; padding: 10px; background-color: white; border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);">
                            <img src="data:image/svg+xml;base64,{image_data}" style="max-width: 100%; height: auto;" />
                        </div>
                        """, unsafe_allow_html=True)
                    except Exception as e:
                        st.error(f"Could not generate visualization: {str(e)}")
            
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
                    
                    # Show the diagram
                    with st.spinner("Generating construct visualization..."):
                        try:
                            # Pass the insert region to highlight it
                            image_data = visualize_construct(construct_result['construct_record'], construct_result.get('insert_region'))
                            
                            # Use HTML component instead of st.image for better SVG rendering
                            st.markdown(f"""
                            <div style="text-align: center; margin: 10px; padding: 10px; background-color: white; border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.1);">
                                <img src="data:image/svg+xml;base64,{image_data}" style="max-width: 100%; height: auto;" />
                            </div>
                            """, unsafe_allow_html=True)
                        except Exception as e:
                            st.error(f"Could not generate visualization: {str(e)}")
                
                # Tab 2: DNA - show the DNA sequence in different formats
                with construct_tabs[1]:
                    sequence_tabs = st.tabs(["FASTA", "GenBank", "Raw Sequence"])
                    
                    with sequence_tabs[0]:
                        st.text(construct_result["fasta_content"])
                        st.download_button(
                            "Download FASTA",
                            construct_result["fasta_content"],
                            file_name=f"{construct_info['gene_name']}_{construct_info['vector_name']}.fasta",
                            mime="text/plain",
                            key="download_fasta_tab"
                        )
                        
                    with sequence_tabs[1]:
                        st.text(construct_result["genbank_content"])
                        st.download_button(
                            "Download GenBank",
                            construct_result["genbank_content"],
                            file_name=f"{construct_info['gene_name']}_{construct_info['vector_name']}.gb",
                            mime="text/plain",
                            key="download_genbank_tab"
                        )
                        
                    with sequence_tabs[2]:
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
                try:
                    with st.sidebar.status("Analyzing construct...", expanded=False) as status:
                        st.sidebar.caption("Examining construct properties...")
                        
                        # Get analysis from the GPT agent
                        analysis = analyze_construct(construct_result['construct_record'], construct_info)
                        
                        # Update status
                        status.update(label="Analysis complete!", state="complete", expanded=True)
                        
                        # Show analysis
                        if "potential_issues" in analysis and analysis["potential_issues"]:
                            st.sidebar.warning("Potential Issues:")
                            for issue in analysis["potential_issues"]:
                                st.sidebar.markdown(f"- {issue}")
                        
                        if "recommendations" in analysis and analysis["recommendations"]:
                            st.sidebar.info("Recommendations:")
                            for rec in analysis["recommendations"]:
                                st.sidebar.markdown(f"- {rec}")
                        
                        if "expression_prediction" in analysis:
                            st.sidebar.markdown(f"**Expression Prediction:** {analysis['expression_prediction']}")
                        
                        if "special_considerations" in analysis and analysis["special_considerations"]:
                            with st.sidebar.expander("Special Considerations"):
                                for consideration in analysis["special_considerations"]:
                                    st.markdown(f"- {consideration}")
                        
                except Exception as e:
                    st.sidebar.error(f"Error analyzing construct: {str(e)}")
                    st.sidebar.info("Analysis not available. Try regenerating the construct.")
            
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