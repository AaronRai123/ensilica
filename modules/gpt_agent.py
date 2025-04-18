"""
LLM Agent module for processing natural language inputs and generating protocols
Using DeepSeek API with hardcoded API key
"""
import json
import os
import requests
from typing import Dict, Any, List, Optional
import logging
from dotenv import load_dotenv
from Bio.SeqRecord import SeqRecord

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load environment variables for API key
load_dotenv()

# Get DeepSeek API key
DEEPSEEK_API_KEY = "sk-0c29dad4f28a4f66ae72e2133c158781"
logger.info("DeepSeek API key found.")

class DeepSeekClient:
    """Client for interacting with the DeepSeek API"""
    
    def __init__(self):
        """Initialize the DeepSeek client with API key"""
        self.api_key = DEEPSEEK_API_KEY
        self.base_url = "https://api.deepseek.com/v1"
        self.headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}"
        }
        
    def chat_completion(self, messages, model="deepseek-chat", temperature=0.7, max_tokens=1000) -> Dict[str, Any]:
        """
        Send a chat completion request to the DeepSeek API
        
        Args:
            messages: List of message objects with role and content
            model: Model name to use
            temperature: Temperature parameter for randomness
            max_tokens: Maximum number of tokens to generate
            
        Returns:
            Response from the API
        """
        if not self.api_key:
            error_msg = "DeepSeek API key not found. Set it as DEEPSEEK_API_KEY environment variable."
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        try:
            logger.info(f"Making DeepSeek API call with model {model}, temp={temperature}")
            url = f"{self.base_url}/chat/completions"
            payload = {
                "model": model,
                "messages": messages,
                "temperature": temperature,
                "max_tokens": max_tokens
            }
            
            response = requests.post(url, headers=self.headers, json=payload)
            response.raise_for_status()  # Raise exception for HTTP errors
            
            result = response.json()
            logger.info(f"DeepSeek API call successful, received {len(result['choices'][0]['message']['content'])} characters")
            return result["choices"][0]["message"]
        except Exception as e:
            logger.error(f"Error calling DeepSeek API: {str(e)}")
            raise ValueError(f"DeepSeek API request failed: {str(e)}")

# Initialize the DeepSeek client
deepseek_client = DeepSeekClient()

def test_api_connection():
    """Test the connection to the DeepSeek API"""
    client = DeepSeekClient()
    try:
        response = client.chat_completion(
            messages=[{"role": "user", "content": "Hello, are you operational?"}],
            max_tokens=10
        )
        return True
    except Exception as e:
        logger.warning(f"DeepSeek API test failed: {str(e)}")
        return False

def process_natural_language(natural_language_input: str) -> Dict[str, Any]:
    """
    Process natural language input to extract structured genetic construct information
    
    Args:
        natural_language_input: Natural language description of genetic construct
        
    Returns:
        Dictionary with structured construct information
    """
    try:
        # Create a more useful prompt with examples
        prompt = f"""
Extract structured information from the following natural language description of a genetic construct.
Return a JSON object with the following fields:
- gene_name: The name of the gene to be cloned (e.g., GFP, mCherry, etc.)
- vector_name: The name of the vector to use (e.g., pET28a, pUC19, etc.)
- promoter_name: The name of the promoter to use (e.g., T7, lac, etc.)
- cloning_method: The cloning method to use (e.g., Gibson Assembly, Restriction Enzyme Cloning)
- insertion_site: Where to insert the gene (e.g., Multiple Cloning Site, specific restriction sites)
- additional_features: Array of additional features to include

EXAMPLES:
Input: "I want to clone GFP into pET28a using Gibson Assembly with a T7 promoter."
Output: {{
  "gene_name": "GFP",
  "vector_name": "pET28a",
  "promoter_name": "T7",
  "cloning_method": "Gibson Assembly",
  "insertion_site": "Multiple Cloning Site",
  "additional_features": []
}}

Input: "Express mCherry in pUC19 using restriction enzyme cloning with EcoRI and HindIII. Use a lac promoter."
Output: {{
  "gene_name": "mCherry",
  "vector_name": "pUC19",
  "promoter_name": "lac",
  "cloning_method": "Restriction Enzyme Cloning",
  "insertion_site": "EcoRI-HindIII",
  "additional_features": []
}}

Now process this input: "{natural_language_input}"
Only respond with the JSON object, nothing else.
"""
        
        client = DeepSeekClient()
        
        # Try using the DeepSeek API
        try:
            response = client.chat_completion(
                messages=[{"role": "user", "content": prompt}],
                temperature=0.1,
                max_tokens=1000
            )
            
            if response and "content" in response:
                # Extract the JSON from the response
                content = response["content"]
                # Find JSON in the response (it might be wrapped in markdown or other text)
                import re
                json_match = re.search(r'{.*}', content, re.DOTALL)
                if json_match:
                    json_str = json_match.group(0)
                    import json
                    result = json.loads(json_str)
                    return result
        except Exception as e:
            logger.warning(f"Error using DeepSeek API: {str(e)}")
            logger.info("Using rule-based fallback for natural language processing")
            
        # Fallback to rule-based parsing if API fails
        return _fallback_natural_language_parser(natural_language_input)
        
    except Exception as e:
        logger.error(f"Error processing natural language: {str(e)}")
        # Return default structure
        return _fallback_natural_language_parser(natural_language_input)

def _fallback_natural_language_parser(text: str) -> Dict[str, Any]:
    """
    Fallback rule-based parser for when the API is not available
    
    Args:
        text: Natural language input
        
    Returns:
        Dictionary with structured construct information
    """
    text = text.lower()
    result = {
        "gene_name": None,
        "vector_name": None,
        "promoter_name": None,
        "cloning_method": "Gibson Assembly",  # Default
        "insertion_site": "Multiple Cloning Site",  # Default
        "additional_features": []
    }
    
    # Look for common genes
    for gene in ["gfp", "mcherry", "rfp", "yfp", "tfeb", "lacz"]:
        if gene in text:
            result["gene_name"] = gene.upper() if gene != "mcherry" else "mCherry"
            break
    
    # Look for common vectors
    for vector in ["pet28a", "puc19", "pbluescript", "pgex", "pcr2.1"]:
        if vector in text:
            # Format vector names correctly
            if vector == "pet28a":
                result["vector_name"] = "pET28a"
            elif vector == "puc19":
                result["vector_name"] = "pUC19"
            elif vector == "pbluescript":
                result["vector_name"] = "pBluescript"
            elif vector == "pgex":
                result["vector_name"] = "pGEX"
            elif vector == "pcr2.1":
                result["vector_name"] = "pCR2.1"
            break
    
    # Look for common promoters
    for promoter in ["t7", "lac", "tac", "tet", "cmv"]:
        if promoter in text:
            result["promoter_name"] = promoter if promoter != "t7" else "T7"
            break
    
    # Look for cloning methods
    if any(method in text for method in ["restriction", "digest", "ecori", "hindiii", "bamhi"]):
        result["cloning_method"] = "Restriction Enzyme Cloning"
        
        # Look for specific restriction sites
        sites = []
        for site in ["ecori", "hindiii", "bamhi", "xhoi", "ncoi"]:
            if site in text:
                sites.append(site.upper())
        
        if sites:
            result["insertion_site"] = "-".join(sites)
    
    # If no gene found, default to GFP
    if result["gene_name"] is None:
        result["gene_name"] = "GFP"
    
    # If no vector found, default to pET28a
    if result["vector_name"] is None:
        result["vector_name"] = "pET28a"
    
    # If no promoter found, default to T7
    if result["promoter_name"] is None:
        result["promoter_name"] = "T7"
    
    return result

def generate_protocol(construct_info: Dict[str, Any], construct_record: SeqRecord = None) -> str:
    """
    Generate a laboratory protocol for the given genetic construct
    
    Args:
        construct_info: Dictionary with construct information
        construct_record: Optional SeqRecord of the construct
        
    Returns:
        String with protocol steps
    """
    try:
        logger.info(f"Generating protocol for {construct_info.get('gene_name', 'unknown')} in {construct_info.get('vector_name', 'unknown')} using {construct_info.get('cloning_method', 'unknown')}")
        
        # Prepare a detailed prompt for protocol generation
        gene_name = construct_info.get("gene_name", "Unknown Gene")
        vector_name = construct_info.get("vector_name", "Unknown Vector")
        cloning_method = construct_info.get("cloning_method", "Unknown Method")
        
        # Create a more specific prompt based on the cloning method
        if "Gibson" in cloning_method:
            method_text = "Gibson Assembly"
            specific_steps = """
1. Design primers with 20-40bp overlaps.
2. PCR amplify the gene of interest with Gibson primers.
3. Linearize the vector with appropriate restriction enzymes or PCR.
4. Set up Gibson Assembly reaction with linearized vector and PCR product.
5. Incubate at 50°C for 15-60 minutes.
6. Transform assembly into competent E. coli.
7. Select transformants on appropriate antibiotic plates.
8. Verify constructs by colony PCR and sequencing.
"""
        elif "Restriction" in cloning_method:
            method_text = "Restriction Enzyme Cloning"
            specific_steps = f"""
1. Digest the vector and gene insert with appropriate restriction enzymes.
2. Dephosphorylate the vector to prevent self-ligation.
3. Gel purify the digested DNA fragments.
4. Set up ligation reaction with T4 DNA Ligase.
5. Incubate ligation at 16°C overnight or room temperature for 1 hour.
6. Transform ligation into competent E. coli.
7. Select transformants on appropriate antibiotic plates.
8. Verify constructs by colony PCR and sequencing.
"""
        else:
            method_text = "Standard Cloning"
            specific_steps = """
1. Prepare the gene insert and vector.
2. Combine the insert and vector using an appropriate cloning method.
3. Transform the resulting construct into competent E. coli.
4. Select transformants on appropriate antibiotic plates.
5. Verify constructs by colony PCR and sequencing.
"""
            
        prompt = f"""
Generate a detailed laboratory protocol for cloning the gene {gene_name} into the vector {vector_name} using {method_text}.
Include specific steps, reagents, temperatures, and incubation times.
The protocol should include:
- Materials and reagents needed
- Step-by-step instructions
- Critical notes and troubleshooting tips

General steps for {method_text}:
{specific_steps}

Return a formatted protocol document suitable for laboratory use.
"""
        
        client = DeepSeekClient()
        
        # Try using the DeepSeek API
        try:
            logger.info("Sending protocol generation request to DeepSeek API")
            response = client.chat_completion(
                messages=[{"role": "user", "content": prompt}],
                temperature=0.3,
                max_tokens=2000
            )
            
            if response and "content" in response:
                logger.info("Successfully generated protocol using DeepSeek API")
                return response["content"]
        except Exception as e:
            logger.error(f"Error generating protocol: {str(e)}")
            logger.info("Using fallback protocol template")
        
        # Fallback to a template if API fails
        logger.info("Using fallback protocol template")
        return _generate_fallback_protocol(construct_info)
    
    except Exception as e:
        logger.error(f"Error generating protocol: {str(e)}")
        return _generate_fallback_protocol(construct_info)

def _generate_fallback_protocol(construct_info: Dict[str, Any]) -> str:
    """Generate a fallback protocol template when the API is unavailable"""
    gene_name = construct_info.get("gene_name", "Unknown Gene")
    vector_name = construct_info.get("vector_name", "Unknown Vector")
    cloning_method = construct_info.get("cloning_method", "Unknown Method")
    promoter_name = construct_info.get("promoter_name", "Unknown Promoter")
    
    if "Gibson" in cloning_method:
        return f"""
# Protocol for Cloning {gene_name} into {vector_name} using Gibson Assembly

## Materials and Reagents
- {gene_name} template DNA
- {vector_name} plasmid DNA
- Q5 High-Fidelity DNA Polymerase (NEB)
- Gibson Assembly Master Mix (NEB)
- Primers for gene amplification with overlaps
- DpnI restriction enzyme
- Competent E. coli cells
- LB medium and agar plates with appropriate antibiotics
- DNA purification kit

## Procedure

### Day 1: PCR and Gibson Assembly

1. Design Gibson Assembly primers with 20-40bp overlaps.
   - Forward primer: [Design specific to {gene_name} with overlap to {vector_name}]
   - Reverse primer: [Design specific to {gene_name} with overlap to {vector_name}]

2. PCR amplify {gene_name} with Gibson primers.
   - PCR reaction (50 μL):
     * 10 μL 5X Q5 Reaction Buffer
     * 1 μL 10 mM dNTPs
     * 2.5 μL Forward Primer (10 μM)
     * 2.5 μL Reverse Primer (10 μM)
     * 1 μL Template DNA (1-10 ng)
     * 0.5 μL Q5 DNA Polymerase
     * 32.5 μL Nuclease-free water
   - PCR cycling:
     * 98°C for 30 seconds
     * 25 cycles of:
       - 98°C for 10 seconds
       - 60°C for 20 seconds
       - 72°C for 30 seconds/kb
     * 72°C for 2 minutes
     * 4°C hold

3. Linearize {vector_name} by PCR or restriction digestion.

4. Purify PCR products using DNA purification kit.

5. Set up Gibson Assembly reaction:
   - 5 μL Gibson Assembly Master Mix (2X)
   - X μL linearized vector (50-100 ng)
   - X μL PCR product (2-3 fold molar excess over vector)
   - Nuclease-free water to 10 μL
   - Incubate at 50°C for 60 minutes

6. Transform 2 μL of the assembly into competent E. coli cells:
   - Add 2 μL assembly to 50 μL competent cells
   - Incubate on ice for 30 minutes
   - Heat shock at 42°C for 30 seconds
   - Return to ice for 2 minutes
   - Add 950 μL SOC medium
   - Incubate at 37°C for 1 hour with shaking
   - Plate 100 μL on LB agar with appropriate antibiotics
   - Incubate plates overnight at 37°C

### Day 2: Colony PCR and Culture

7. Pick 5-10 colonies for verification by colony PCR.
8. Inoculate positive colonies into LB medium with antibiotics.
9. Grow overnight at 37°C with shaking.

### Day 3: Plasmid Extraction and Verification

10. Extract plasmid DNA using a miniprep kit.
11. Verify the construct by Sanger sequencing.

## Notes
- For optimal Gibson Assembly, the overlapping regions should be 20-40 bp with similar melting temperatures.
- The {promoter_name} promoter should drive expression of {gene_name}.
- Verify the final construct by both restriction digestion and sequencing.
"""
    else:  # Restriction Enzyme Cloning
        return f"""
# Protocol for Cloning {gene_name} into {vector_name} using Restriction Enzyme Cloning

## Materials and Reagents
- {gene_name} template DNA
- {vector_name} plasmid DNA
- Restriction enzymes (based on your insertion site)
- T4 DNA Ligase and 10X buffer
- Antarctic Phosphatase
- Q5 High-Fidelity DNA Polymerase
- Gel extraction kit
- Competent E. coli cells
- LB medium and agar plates with appropriate antibiotics

## Procedure

### Day 1: Restriction Digestion and Ligation

1. PCR amplify {gene_name} with primers containing restriction sites.
   - PCR reaction (50 μL):
     * 10 μL 5X Q5 Reaction Buffer
     * 1 μL 10 mM dNTPs
     * 2.5 μL Forward Primer (10 μM)
     * 2.5 μL Reverse Primer (10 μM)
     * 1 μL Template DNA (1-10 ng)
     * 0.5 μL Q5 DNA Polymerase
     * 32.5 μL Nuclease-free water
   - PCR cycling:
     * 98°C for 30 seconds
     * 25 cycles of:
       - 98°C for 10 seconds
       - 60°C for 20 seconds
       - 72°C for 30 seconds/kb
     * 72°C for 2 minutes
     * 4°C hold

2. Digest PCR product and vector with appropriate restriction enzymes.
   - Reaction (50 μL):
     * 5 μL 10X Restriction Buffer
     * 1-5 μg DNA
     * 1 μL each Restriction Enzyme
     * Nuclease-free water to 50 μL
   - Incubate at 37°C for 1-2 hours

3. Dephosphorylate the vector.
   - Add 1 μL Antarctic Phosphatase and 5 μL 10X buffer to the digested vector
   - Incubate at 37°C for 30 minutes
   - Heat inactivate at 80°C for 2 minutes

4. Gel purify digested DNA fragments.

5. Set up ligation reaction:
   - 2 μL 10X T4 DNA Ligase Buffer
   - 50 ng digested vector
   - Insert DNA (3:1 molar ratio insert:vector)
   - 1 μL T4 DNA Ligase
   - Nuclease-free water to 20 μL
   - Incubate at 16°C overnight or room temperature for 1 hour

6. Transform 5 μL of the ligation into competent E. coli cells:
   - Add 5 μL ligation to 50 μL competent cells
   - Incubate on ice for 30 minutes
   - Heat shock at 42°C for 30 seconds
   - Return to ice for 2 minutes
   - Add 950 μL SOC medium
   - Incubate at 37°C for 1 hour with shaking
   - Plate 100 μL on LB agar with appropriate antibiotics
   - Incubate plates overnight at 37°C

### Day 2: Colony PCR and Culture

7. Pick 5-10 colonies for verification by colony PCR.
8. Inoculate positive colonies into LB medium with antibiotics.
9. Grow overnight at 37°C with shaking.

### Day 3: Plasmid Extraction and Verification

10. Extract plasmid DNA using a miniprep kit.
11. Verify the construct by restriction digestion and Sanger sequencing.

## Notes
- The {promoter_name} promoter should drive expression of {gene_name}.
- Ensure insert is in the correct orientation after ligation.
- Include appropriate controls (vector-only ligation, etc.)
"""

def analyze_construct(construct_info, construct_record=None):
    """
    Analyze a genetic construct and provide feedback
    
    Args:
        construct_info: Dictionary with construct information or SeqRecord object
        construct_record: Optional SeqRecord of the construct
        
    Returns:
        Dictionary with analysis results
    """
    try:
        # Handle the case where construct_info is a SeqRecord
        if hasattr(construct_info, 'seq') and not hasattr(construct_info, 'get'):
            if construct_record is None:
                construct_record = construct_info
            # Create a default construct_info dictionary
            construct_info = {
                'gene_name': getattr(construct_record, 'id', 'Unknown Gene').split('_')[0] if '_' in getattr(construct_record, 'id', '') else getattr(construct_record, 'id', 'Unknown Gene'),
                'vector_name': getattr(construct_record, 'id', 'Unknown Vector').split('_')[1] if '_' in getattr(construct_record, 'id', '') and len(getattr(construct_record, 'id', '').split('_')) > 1 else 'Unknown Vector',
                'promoter_name': 'Unknown Promoter',
                'cloning_method': 'Unknown Method'
            }
            
            # Try to extract more info from features if available
            if hasattr(construct_record, 'features'):
                for feature in construct_record.features:
                    if hasattr(feature, 'type') and feature.type == 'promoter' and hasattr(feature, 'qualifiers'):
                        qualifiers = feature.qualifiers
                        if isinstance(qualifiers, dict) and 'label' in qualifiers:
                            labels = qualifiers['label']
                            if isinstance(labels, list) and labels:
                                construct_info['promoter_name'] = labels[0]
        
        logger.info(f"Analyzing construct {construct_info.get('gene_name', 'unknown')} in {construct_info.get('vector_name', 'unknown')}")
        
        # Prepare a prompt for analysis
        gene_name = construct_info.get("gene_name", "Unknown Gene")
        vector_name = construct_info.get("vector_name", "Unknown Vector")
        promoter_name = construct_info.get("promoter_name", "Unknown Promoter")
        cloning_method = construct_info.get("cloning_method", "Unknown Method")
        
        # Additional information if available
        seq_length = "Unknown"
        has_origin = False
        has_promoter = False
        
        # Safely extract information from construct_record if available
        if construct_record is not None:
            logger.info("Extracting information from construct record")
            if hasattr(construct_record, 'seq'):
                seq_length = len(construct_record.seq)
                logger.info(f"Construct sequence length: {seq_length} bp")
            
            if hasattr(construct_record, 'features'):
                logger.info(f"Construct has {len(construct_record.features)} features")
                for feature in construct_record.features:
                    # Check for origin of replication
                    if hasattr(feature, 'type'):
                        logger.info(f"Found feature type: {feature.type}")
                        if 'origin' in feature.type.lower():
                            has_origin = True
                    
                    # Check for promoter
                    if hasattr(feature, 'type'):
                        if 'promoter' in feature.type.lower():
                            has_promoter = True
                    
                    # Also check qualifiers if available
                    if hasattr(feature, 'qualifiers'):
                        qualifiers = feature.qualifiers
                        if isinstance(qualifiers, dict):
                            # Check labels
                            if 'label' in qualifiers:
                                labels = qualifiers['label']
                                if isinstance(labels, list) and labels:
                                    label = labels[0].lower()
                                    logger.info(f"Found feature label: {label}")
                                    if 'ori' in label or 'origin' in label:
                                        has_origin = True
                                    if 'promot' in label:
                                        has_promoter = True
                logger.info(f"Construct has origin: {has_origin}, has promoter: {has_promoter}")
        
        prompt = f"""
Analyze the following genetic construct and provide feedback on its design and potential issues:

Gene: {gene_name}
Vector: {vector_name}
Promoter: {promoter_name}
Cloning Method: {cloning_method}
Construct Length: {seq_length} bp
Contains Origin: {has_origin}
Contains Promoter: {has_promoter}

Provide feedback on:
1. Suitability of the chosen vector for this gene
2. Expected expression levels
3. Potential cloning complications
4. Recommendations for validation
5. Potential applications

Return your analysis as a JSON object with these 5 sections.
"""
        
        client = DeepSeekClient()
        
        # Try using the DeepSeek API
        try:
            logger.info("Sending construct analysis request to DeepSeek API")
            response = client.chat_completion(
                messages=[{"role": "user", "content": prompt}],
                temperature=0.3,
                max_tokens=1000
            )
            
            if response and "content" in response:
                logger.info("Successfully generated analysis using DeepSeek API")
                # Extract the JSON from the response
                content = response["content"]
                import re
                json_match = re.search(r'{.*}', content, re.DOTALL)
                if json_match:
                    json_str = json_match.group(0)
                    import json
                    return json.loads(json_str)
        except Exception as e:
            logger.error(f"Error analyzing construct: {str(e)}")
            logger.info("Using fallback analysis logic")
        
        # Fallback to a template if API fails
        logger.info("Using fallback analysis template")
        return _generate_fallback_analysis(construct_info, construct_record)
    
    except Exception as e:
        logger.error(f"Error analyzing construct: {str(e)}")
        return _generate_fallback_analysis(construct_info, construct_record)

def _generate_fallback_analysis(construct_info, construct_record = None) -> Dict[str, Any]:
    """Generate fallback analysis when the API is unavailable"""
    # Handle the case where construct_info is a SeqRecord
    if hasattr(construct_info, 'seq') and not hasattr(construct_info, 'get'):
        if construct_record is None:
            construct_record = construct_info
        # Create a default construct_info dictionary
        construct_info = {
            'gene_name': getattr(construct_record, 'id', 'Unknown Gene').split('_')[0] if '_' in getattr(construct_record, 'id', '') else getattr(construct_record, 'id', 'Unknown Gene'),
            'vector_name': getattr(construct_record, 'id', 'Unknown Vector').split('_')[1] if '_' in getattr(construct_record, 'id', '') and len(getattr(construct_record, 'id', '').split('_')) > 1 else 'Unknown Vector',
            'promoter_name': 'Unknown Promoter',
            'cloning_method': 'Unknown Method'
        }
    
    gene_name = construct_info.get("gene_name", "Unknown Gene")
    vector_name = construct_info.get("vector_name", "Unknown Vector")
    promoter_name = construct_info.get("promoter_name", "Unknown Promoter")
    cloning_method = construct_info.get("cloning_method", "Unknown Method")
    
    # Default analysis
    analysis = {
        "vector_suitability": "The selected vector appears to be appropriate for this construct.",
        "expression_levels": "Expression levels are expected to be moderate to high depending on growth conditions.",
        "cloning_complications": "Standard cloning challenges may apply. Verify the final construct by sequencing.",
        "validation_recommendations": "Verify by restriction digestion and sequencing. Test expression under various conditions.",
        "potential_applications": "This construct could be used for protein expression and purification studies."
    }
    
    # Gene-specific feedback
    if gene_name == "GFP" or gene_name == "mCherry" or gene_name == "RFP":
        analysis["expression_levels"] = f"The {gene_name} fluorescent protein should be visibly detectable if expression is successful."
        analysis["validation_recommendations"] += f" Fluorescence microscopy can be used to verify {gene_name} expression."
        analysis["potential_applications"] = f"This {gene_name} construct can be used for visualization studies, protein localization, or as a reporter gene."
    
    # Vector-specific feedback
    if vector_name == "pET28a":
        analysis["vector_suitability"] = "pET28a is a high-expression T7 promoter-based vector with a His-tag, suitable for bacterial protein expression."
        analysis["expression_levels"] = "High expression levels in E. coli DE3 strains upon IPTG induction."
        analysis["cloning_complications"] = "Ensure the insert is in-frame with the His-tag if protein purification is intended."
    elif vector_name == "pUC19":
        analysis["vector_suitability"] = "pUC19 is a high-copy cloning vector with a lac promoter, suitable for DNA propagation."
        analysis["expression_levels"] = "The lac promoter provides moderate expression levels in E. coli."
        analysis["cloning_complications"] = "The high copy number may be challenging if the insert is toxic to cells."
    elif vector_name == "pBluescript":
        analysis["vector_suitability"] = "pBluescript is a general-purpose cloning vector with blue-white screening capability."
        analysis["expression_levels"] = "Expression levels are typically low without an added strong promoter."
        analysis["cloning_complications"] = "Use blue-white screening with X-gal to identify recombinant colonies."
    
    # Promoter-specific feedback
    if promoter_name == "T7":
        analysis["expression_levels"] = "The T7 promoter provides very high expression levels but requires a DE3 strain for T7 RNA polymerase."
    elif promoter_name == "lac":
        analysis["expression_levels"] = "The lac promoter provides moderate, inducible expression with IPTG."
    
    # Cloning method feedback
    if "Gibson" in cloning_method:
        analysis["cloning_complications"] = "Gibson Assembly is generally efficient but requires careful primer design with appropriate overlaps."
    elif "Restriction" in cloning_method:
        analysis["cloning_complications"] = "Restriction enzyme cloning may require optimization of digestion and ligation conditions. Check for internal restriction sites."
    
    return analysis

# Run a connection test when the module is imported
if __name__ == "__main__":
    test_api_connection() 