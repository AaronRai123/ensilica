#!/bin/bash

# Check if the DeepSeek API key is set
if [ -z "$DEEPSEEK_API_KEY" ]; then
    echo "Warning: DEEPSEEK_API_KEY environment variable is not set."
    echo "The app will run with limited functionality."
    echo "Set it with: export DEEPSEEK_API_KEY=your_api_key_here"
fi

# Run the setup script to ensure directories exist
python setup.py

# Launch the Streamlit app
streamlit run app.py 