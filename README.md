# Ensilica

## AI-Powered Genetic Construct Designer

Ensilica is an AI-powered Copilot that allows bioengineers to design genetic constructs using natural language. Simply describe what you want, and Ensilica will generate the DNA sequence, primers, and a lab-ready protocol.

### 🧪 Features

- **Natural Language Input**: Describe your desired genetic construct in plain English
- **DNA Sequence Generation**: Automatically create DNA constructs by inserting genes into vectors
- **Primer Design**: Generate PCR primers for your construct using Primer3
- **Lab Protocol**: Receive a detailed, step-by-step protocol for your experiment
- **Downloadable Files**: Export everything you need in standard formats (FASTA, GenBank, etc.)

### 📦 Requirements

- Python 3.8+
- Streamlit
- Biopython
- Pandas
- OpenAI API key (or DeepSeek API key)
- primer3-py (optional)

### 🛠️ Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/ensilica.git
   cd ensilica
   ```

2. Install dependencies:
   ```
   pip install -r requirements.txt
   ```

3. Set up your API key:
   ```
   export OPENAI_API_KEY=your_api_key_here
   ```
   Or for DeepSeek:
   ```
   export DEEPSEEK_API_KEY=your_api_key_here
   ```

### 🚀 Usage

1. Start the Streamlit app:
   ```
   streamlit run app.py
   ```

2. Open your browser at `http://localhost:8501`

3. Enter a natural language description of your genetic construct, such as:
   ```
   Insert mCherry into pUC19 