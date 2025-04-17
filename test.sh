#!/bin/bash

# Colors for terminal output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Ensilica Testing Script${NC}"
echo -e "${YELLOW}====================${NC}\n"

# Check if DeepSeek API key is set
if [ -z "$DEEPSEEK_API_KEY" ]; then
    echo -e "${RED}Warning: DEEPSEEK_API_KEY environment variable is not set.${NC}"
    echo -e "API-dependent tests will fail. Set the key with:

export DEEPSEEK_API_KEY=your_api_key_here\n"
fi 