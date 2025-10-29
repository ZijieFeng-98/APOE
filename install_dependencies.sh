#!/bin/bash
# Install Python dependencies for clinical analysis
# Run this from within WSL

echo "=================================="
echo "Installing Clinical Analysis Dependencies"
echo "=================================="

# Update package lists
echo "📦 Updating package lists..."
sudo apt update

# Install pip if not already installed
echo "📦 Installing pip..."
sudo apt install -y python3-pip

# Install required Python packages
echo "📦 Installing Python packages..."
pip3 install --user pandas matplotlib seaborn lifelines openpyxl

# Verify installation
echo ""
echo "✅ Verifying installation..."
python3 -c "import pandas; print(f'✓ pandas {pandas.__version__}')"
python3 -c "import matplotlib; print(f'✓ matplotlib {matplotlib.__version__}')"
python3 -c "import seaborn; print(f'✓ seaborn {seaborn.__version__}')"
python3 -c "import lifelines; print(f'✓ lifelines {lifelines.__version__}')"
python3 -c "import openpyxl; print(f'✓ openpyxl {openpyxl.__version__}')"

echo ""
echo "=================================="
echo "✅ Installation complete!"
echo "=================================="
echo ""
echo "Now run: python3 run_clinical_analysis.py"

