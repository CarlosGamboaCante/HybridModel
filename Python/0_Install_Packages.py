# ============================================================
# HybridModel Tool
# Install required Python packages
#
# Compatible with:
#   Python 3.11
#
# The package versions below are recommended to ensure
# compatibility with the current version of HybridModel.
# ============================================================

import sys
import subprocess

print("Installing required packages for HybridModel...\n")

packages = [
    "numpy==1.26.4",
    "pandas==1.5.3",
    "scipy==1.11.4",
    "matplotlib==3.8.4",
]

for package in packages:
    print(f"Installing {package}...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

print("\n==========================================")
print("All packages were installed successfully.")
print("HybridModel is ready to use.")
print("==========================================")
