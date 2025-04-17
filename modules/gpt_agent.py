"""
LLM Agent module for processing natural language inputs and generating protocols
Using DeepSeek API with hardcoded API key
"""
import json
import os
import requests
from typing import Dict, Any, List, Optional
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Hardcoded DeepSeek API key - this will actually call the DeepSeek API
# Replace this with your actual DeepSeek API key
DEEPSEEK_API_KEY = "sk-31d8b2d0bd99415294517c3b7afb71cf"
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"

class DeepSeekClient:
    """Client for interacting with DeepSeek API"""
    
    def __init__(self):
        self.api_key = DEEPSEEK_API_KEY
        self.headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}"
        }
    
    def chat_completion(self, 
                        messages: List[Dict[str, str]], 
                        model: str = "deepseek-chat",
                        temperature: float = 0.7,
                        max_tokens: int = 2000) -> Dict[str, Any]:
        """
        Send a chat completion request to DeepSeek API
        
        Args:
            messages: List of message dictionaries with 'role' and 'content' keys
            model: DeepSeek model to use
            temperature: Sampling temperature (0.0 to 1.0)
            max_tokens: Maximum number of tokens to generate
            
        Returns:
            DeepSeek API response
        """
        payload = {
            "model": model,
            "messages": messages,
            "temperature": temperature,
            "max_tokens": max_tokens
        }
        
        try:
            response = requests.post(
                DEEPSEEK_API_URL,
                headers=self.headers,
                json=payload
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"DeepSeek API request failed: {str(e)}")
            if response := getattr(e, 'response', None):
                logger.error(f"Response status: {response.status_code}")
                logger.error(f"Response body: {response.text}")
            raise

# Initialize DeepSeek client
deepseek_client = DeepSeekClient()

def process_natural_language(input_text: str) -> Dict[str, Any]:
    """
    Process natural language input to extract structured information about the genetic construct.
    
    Args:
        input_text: Natural language description of the desired genetic construct
        
    Returns:
        Dictionary containing structured information about the construct
    """
    # Prepare the prompt for DeepSeek
    prompt = f"""
    Extract structured information from the following description of a genetic construct:
    "{input_text}"
    
    Return a JSON object with the following fields:
    - gene_name: The name of the gene to be inserted (e.g., mCherry, GFP)
    - vector_name: The name of the vector backbone (e.g., pUC19, pET28a)
    - promoter_name: The name of the promoter (e.g., lac, T7)
    - cloning_method: The method to be used for cloning (e.g., Gibson Assembly, restriction digestion)
    - insertion_site: Where the gene should be inserted, if specified
    - additional_features: Any other features mentioned (e.g., tags, selection markers)
    
    If any field is not specified, use null or provide a sensible default.
    Format the response as valid JSON only, with no additional text.
    """
    
    try:
        # Call the DeepSeek API
        response = deepseek_client.chat_completion(
            messages=[
                {"role": "system", "content": "You are a biological engineering assistant that extracts structured information from text and returns valid JSON only."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.2  # Lower temperature for more deterministic outputs
        )
        
        # Parse the response as JSON
        content = response.get("choices", [{}])[0].get("message", {}).get("content", "{}")
        result = json.loads(content)
        return result
    
    except Exception as e:
        # Log the error
        logger.error(f"Error processing natural language: {str(e)}")
        
        # Fallback to basic parsing if API call fails
        logger.info("Using fallback parsing logic")
        
        # Basic parsing logic
        result = {
            "gene_name": "GFP" if "GFP" in input_text else ("mCherry" if "mCherry" in input_text else "unknown_gene"),
            "vector_name": "pUC19" if "pUC19" in input_text else ("pET28a" if "pET28a" in input_text else "pUC19"),
            "promoter_name": "lac" if "lac" in input_text else ("T7" if "T7" in input_text else "lac"),
            "cloning_method": "Gibson Assembly" if "Gibson" in input_text else "Restriction Enzyme Cloning",
            "insertion_site": "Multiple Cloning Site",
            "additional_features": []  # Empty list instead of None
        }
        
        return result

def generate_protocol(construct_info: Dict[str, Any], primer_info: Dict[str, Any]) -> str:
    """
    Generate a laboratory protocol for the genetic construct based on the extracted information.
    
    Args:
        construct_info: Dictionary containing structured information about the construct
        primer_info: Dictionary containing primer information
        
    Returns:
        Markdown-formatted laboratory protocol
    """
    # Get key information
    gene_name = construct_info.get("gene_name", "unknown_gene")
    vector_name = construct_info.get("vector_name", "unknown_vector")
    promoter_name = construct_info.get("promoter_name", "unknown_promoter")
    cloning_method = construct_info.get("cloning_method", "Gibson Assembly")
    
    # Format primers for the protocol
    forward_primer = primer_info['primers'][0]['sequence']
    reverse_primer = primer_info['primers'][1]['sequence']
    
    # Prepare prompt for DeepSeek
    prompt = f"""
    Generate a detailed laboratory protocol for inserting {gene_name} into {vector_name} under the {promoter_name} promoter using {cloning_method}.
    
    Include the following information:
    - Materials needed
    - Step-by-step procedure
    - PCR conditions
    - Assembly conditions
    - Transformation protocol
    - Verification steps
    
    Use these PCR primers:
    Forward: {forward_primer}
    Reverse: {reverse_primer}
    
    Format the output as Markdown.
    """
    
    try:
        # Call the DeepSeek API
        response = deepseek_client.chat_completion(
            messages=[
                {"role": "system", "content": "You are a biological engineering assistant that generates detailed, accurate laboratory protocols."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7,  # More creativity for protocol generation
            max_tokens=2500   # Protocols can be lengthy
        )
        
        # Extract the protocol text
        protocol = response.get("choices", [{}])[0].get("message", {}).get("content", "")
        return protocol
    
    except Exception as e:
        # Log the error
        logger.error(f"Error generating protocol: {str(e)}")
        
        # Fallback to template protocol if API call fails
        logger.info("Using fallback protocol template")
        
        # Generate a template protocol based on cloning method
        if cloning_method == "Gibson Assembly":
            protocol = f"""
# Laboratory Protocol: {gene_name} insertion into {vector_name} under {promoter_name} promoter

## Materials
- {vector_name} plasmid DNA (100 ng/μL)
- {gene_name} template DNA (100 ng/μL)
- Forward primer: {forward_primer}
- Reverse primer: {reverse_primer}
- Q5 High-Fidelity DNA Polymerase
- Gibson Assembly Master Mix
- Competent E. coli cells
- LB agar plates with appropriate antibiotic

## Procedure

### 1. PCR Amplification
1. Set up PCR reaction to amplify {gene_name} with added overlaps for Gibson Assembly
2. Purify PCR product

### 2. Vector Preparation
1. Linearize {vector_name} by PCR or restriction digestion
2. Purify linearized vector

### 3. Gibson Assembly
1. Mix insert and vector in a 3:1 molar ratio with Gibson Assembly Master Mix
2. Incubate at 50°C for 60 minutes

### 4. Transformation
1. Transform into competent cells
2. Plate on selective media
3. Incubate overnight

### 5. Verification
1. Screen colonies by PCR
2. Verify by sequencing

For detailed steps, please refer to NEB's Gibson Assembly protocol.
            """
        else:
            protocol = f"""
# Laboratory Protocol: {gene_name} insertion into {vector_name}

## Materials
- {vector_name} plasmid DNA
- {gene_name} DNA
- Forward primer: {forward_primer}
- Reverse primer: {reverse_primer}
- Restriction enzymes
- T4 DNA Ligase
- Competent E. coli cells
- LB agar plates with appropriate antibiotic

## Procedure

### 1. Digest Vector and Insert
1. Digest {vector_name} with appropriate restriction enzymes
2. Digest {gene_name} with the same restriction enzymes
3. Purify digested products

### 2. Ligation
1. Mix vector and insert in a 1:3 molar ratio with T4 DNA Ligase
2. Incubate at 16°C overnight

### 3. Transformation
1. Transform into competent cells
2. Plate on selective media
3. Incubate overnight

### 4. Verification
1. Screen colonies by PCR
2. Verify by restriction digestion
3. Confirm by sequencing

For detailed steps, please refer to standard molecular cloning protocols.
            """
        
        return protocol

def test_api_connection() -> bool:
    """Test the connection to the DeepSeek API"""
    if not DEEPSEEK_API_KEY:
        logger.error("No DeepSeek API key provided")
        return False
    
    try:
        # Simple test query
        response = deepseek_client.chat_completion(
            messages=[
                {"role": "user", "content": "Hello, are you working?"}
            ],
            max_tokens=10
        )
        logger.info("DeepSeek API connection successful")
        return True
    except Exception as e:
        logger.error(f"DeepSeek API connection failed: {str(e)}")
        return False

def analyze_construct(construct_record, construct_info):
    """
    Analyze a genetic construct and provide feedback.
    
    Args:
        construct_record: SeqRecord object containing the construct
        construct_info: Dictionary with construct specifications
        
    Returns:
        Dictionary with analysis results including potential issues,
        recommendations, expression prediction, and special considerations
    """
    # Extract basic information
    gene_name = construct_info.get("gene_name", "unknown")
    vector_name = construct_info.get("vector_name", "unknown")
    promoter_name = construct_info.get("promoter_name", "unknown")
    cloning_method = construct_info.get("cloning_method", "unknown")
    
    # Get some basic sequence properties
    seq_length = len(construct_record.seq) if construct_record and hasattr(construct_record, 'seq') else 0
    feature_count = len(construct_record.features) if construct_record and hasattr(construct_record, 'features') else 0
    
    # Prepare prompt for DeepSeek
    prompt = f"""
    Analyze the following genetic construct and provide feedback:
    
    Construct details:
    - Gene: {gene_name}
    - Vector: {vector_name}
    - Promoter: {promoter_name}
    - Cloning method: {cloning_method}
    - Total length: {seq_length} bp
    - Number of features: {feature_count}
    
    Provide analysis in JSON format with these fields:
    - potential_issues: List of any potential issues with the construct
    - recommendations: List of recommendations for optimization
    - expression_prediction: Estimated expression level (High/Medium/Low)
    - special_considerations: Any special considerations for this specific construct
    
    Format as valid JSON only.
    """
    
    try:
        # Call the DeepSeek API
        response = deepseek_client.chat_completion(
            messages=[
                {"role": "system", "content": "You are a biological engineering assistant that analyzes genetic constructs and provides expert feedback."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.5,
            max_tokens=1000
        )
        
        # Extract the analysis
        content = response.get("choices", [{}])[0].get("message", {}).get("content", "{}")
        analysis = json.loads(content)
        return analysis
    
    except Exception as e:
        # Log the error
        logger.error(f"Error analyzing construct: {str(e)}")
        
        # Fallback to basic analysis if API call fails
        logger.info("Using fallback analysis logic")
        
        # Generate basic analysis
        if promoter_name == "T7" and vector_name != "pET28a":
            potential_issues = ["T7 promoter is typically used with pET vectors; may not function properly in " + vector_name]
            exp_level = "Low"
        elif promoter_name == "T7" and gene_name in ["GFP", "mCherry", "RFP"]:
            potential_issues = []
            exp_level = "High"
        else:
            potential_issues = []
            exp_level = "Medium"
            
        # Common recommendations based on gene and vector
        if gene_name in ["GFP", "mCherry", "RFP"]:
            recommendations = ["Consider codon optimization for your expression host", 
                              "Add a strong RBS upstream of your coding sequence"]
        else:
            recommendations = ["Verify insert orientation after cloning", 
                              "Consider sequencing the entire construct before expression"]
                              
        # Special considerations
        if "His tag" in str(construct_info):
            special_considerations = "The His-tag may affect protein folding or function; consider testing both N and C-terminal versions"
        else:
            special_considerations = "No special considerations identified"
            
        return {
            "potential_issues": potential_issues,
            "recommendations": recommendations,
            "expression_prediction": exp_level,
            "special_considerations": special_considerations
        }

# Run a connection test when the module is imported
if __name__ == "__main__":
    test_api_connection() 