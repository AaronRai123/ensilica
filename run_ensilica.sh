#!/bin/bash

# Colors for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=======================================${NC}"
echo -e "${GREEN}  Ensilica - AI Genetic Construct Designer${NC}"
echo -e "${GREEN}=======================================${NC}"

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: Python 3 is not installed. Please install Python 3 and try again.${NC}"
    exit 1
fi

# Check if required packages are installed
if ! pip list | grep -q "streamlit"; then
    echo -e "${YELLOW}Installing required packages...${NC}"
    pip install -r requirements.txt
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}Failed to install required packages. Let's try specific versions...${NC}"
        pip install streamlit==1.24.0 biopython==1.81 pandas==1.5.3 matplotlib==3.7.1 numpy==1.24.3 requests==2.31.0 python-dotenv==1.0.0
    fi
fi

# Check if DeepSeek API key is set
if [ -z "$DEEPSEEK_API_KEY" ]; then
    echo -e "${YELLOW}Warning: DeepSeek API key is not set in environment.${NC}"
    echo -e "${YELLOW}You can set it in the app or set it now:${NC}"
    read -p "Enter DeepSeek API Key (leave empty to skip): " api_key
    
    if [ ! -z "$api_key" ]; then
        export DEEPSEEK_API_KEY="$api_key"
        echo -e "${GREEN}DeepSeek API key set!${NC}"
    else
        echo -e "${YELLOW}No API key provided. The app will use fallback mechanisms.${NC}"
    fi
fi

# Run the app
echo -e "${GREEN}Starting Ensilica...${NC}"
streamlit run app_fixed.py

# Exit gracefully
echo -e "${GREEN}Ensilica session ended.${NC}" 