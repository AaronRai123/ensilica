import os
import shutil

# Create necessary directories
os.makedirs("data", exist_ok=True)
os.makedirs("modules", exist_ok=True)

# Ensure the empty __init__.py file exists in modules
with open("modules/__init__.py", "w") as f:
    f.write("# This file makes the modules directory a Python package")

print("Setup complete. Directory structure created.") 