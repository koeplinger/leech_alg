import sys
import pathlib

# Make src/ importable from any test file under python_project/
sys.path.insert(0, str(pathlib.Path(__file__).parent / "src"))
