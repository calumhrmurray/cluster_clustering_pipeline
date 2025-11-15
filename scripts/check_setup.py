#!/usr/bin/env python3
"""
Check that the pipeline is set up correctly and all dependencies are installed.
"""

import sys
from pathlib import Path

def check_python_version():
    """Check Python version."""
    version = sys.version_info
    print(f"Python version: {version.major}.{version.minor}.{version.micro}")
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        print("  ❌ Python 3.8+ required")
        return False
    print("  ✓ Python version OK")
    return True

def check_dependencies():
    """Check required packages are installed."""
    packages = {
        'astropy': 'astropy',
        'numpy': 'numpy',
        'treecorr': 'treecorr',
        'yaml': 'pyyaml',
        'matplotlib': 'matplotlib (optional, for plotting)'
    }

    all_ok = True
    for module, name in packages.items():
        try:
            __import__(module)
            print(f"  ✓ {name}")
        except ImportError:
            print(f"  ❌ {name} not found")
            all_ok = False

    return all_ok

def check_directory_structure():
    """Check pipeline directory structure."""
    pipeline_dir = Path(__file__).parent.parent

    required_dirs = ['src', 'scripts', 'configs']
    required_files = [
        'src/__init__.py',
        'src/catalogue.py',
        'src/filters.py',
        'src/clustering.py',
        'src/config.py',
        'scripts/run_analysis.py',
        'scripts/generate_jobs.py',
        'README.md'
    ]

    all_ok = True

    # Check directories
    for dir_name in required_dirs:
        dir_path = pipeline_dir / dir_name
        if dir_path.exists():
            print(f"  ✓ {dir_name}/")
        else:
            print(f"  ❌ {dir_name}/ not found")
            all_ok = False

    # Check files
    for file_name in required_files:
        file_path = pipeline_dir / file_name
        if file_path.exists():
            print(f"  ✓ {file_name}")
        else:
            print(f"  ❌ {file_name} not found")
            all_ok = False

    return all_ok

def check_imports():
    """Check pipeline modules can be imported."""
    sys.path.insert(0, str(Path(__file__).parent.parent))

    modules = [
        'src',
        'src.catalogue',
        'src.filters',
        'src.clustering',
        'src.shapes',
        'src.weights',
        'src.config'
    ]

    all_ok = True
    for module in modules:
        try:
            __import__(module)
            print(f"  ✓ {module}")
        except Exception as e:
            print(f"  ❌ {module}: {e}")
            all_ok = False

    return all_ok

def main():
    print("="*70)
    print("CLUSTER CLUSTERING PIPELINE - SETUP CHECK")
    print("="*70)

    checks = [
        ("Python Version", check_python_version),
        ("Dependencies", check_dependencies),
        ("Directory Structure", check_directory_structure),
        ("Module Imports", check_imports)
    ]

    results = {}
    for name, check_func in checks:
        print(f"\n{name}:")
        results[name] = check_func()

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    all_passed = all(results.values())

    for name, passed in results.items():
        status = "✓ PASS" if passed else "❌ FAIL"
        print(f"{name}: {status}")

    if all_passed:
        print("\n✓ All checks passed! Pipeline is ready to use.")
        print("\nNext steps:")
        print("1. Create a config: python scripts/create_config.py")
        print("2. Test the pipeline: python scripts/run_analysis.py --config config.yaml --test")
        print("3. See QUICKSTART.md for more details")
        return 0
    else:
        print("\n❌ Some checks failed. Please install missing dependencies:")
        print("   pip install -r requirements.txt")
        return 1

if __name__ == '__main__':
    sys.exit(main())
