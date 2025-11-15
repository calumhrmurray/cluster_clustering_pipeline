# Installation Guide

## Prerequisites

- Python 3.8 or higher
- Access to SLURM cluster (for batch job submission)
- Cluster and galaxy catalogues in FITS format

## Installation Steps

### 1. Navigate to Pipeline Directory

```bash
cd /sps/euclid/Users/cmurray/euclid_cluster_clustering/cluster_clustering_pipeline
```

### 2. Install Dependencies

**Option A: Using pip**
```bash
pip install -r requirements.txt
```

**Option B: Using conda (recommended)**
```bash
# Create new environment
conda create -n cluster_pipeline python=3.11
conda activate cluster_pipeline

# Install packages
pip install -r requirements.txt

# Or install via conda
conda install -c conda-forge astropy numpy pyyaml matplotlib
pip install treecorr
```

**Option C: On SLURM cluster with modules**
```bash
# Load Python module
module load python/3.11

# Create virtual environment
python -m venv ~/envs/cluster_pipeline
source ~/envs/cluster_pipeline/bin/activate

# Install packages
pip install -r requirements.txt
```

### 3. Verify Installation

```bash
python scripts/check_setup.py
```

You should see all checks pass. If treecorr installation fails, see Troubleshooting below.

### 4. Test Pipeline

```bash
# Create a test config
python scripts/create_config.py --output test_config.yaml

# Edit test_config.yaml with your catalogue paths
# Then run a test
python scripts/run_analysis.py --config test_config.yaml --test
```

## Environment Setup for SLURM

Add to your config file:

```yaml
slurm:
  # Load required modules
  modules:
    - python/3.11
  # Activate conda environment
  conda_env: cluster_pipeline
  # OR use virtualenv
  # venv: /path/to/your/venv
```

The job scripts will automatically load modules and activate environments.

## Directory Structure

After installation, you should have:

```
cluster_clustering_pipeline/
├── src/                     # Core Python modules
│   ├── __init__.py
│   ├── catalogue.py
│   ├── clustering.py
│   ├── config.py
│   ├── filters.py
│   ├── shapes.py
│   └── weights.py
├── scripts/                 # Command-line tools
│   ├── check_setup.py
│   ├── create_config.py
│   ├── generate_jobs.py
│   ├── plot_results.py
│   ├── run_analysis.py
│   └── run_single_job.py
├── configs/                 # Example configurations
│   ├── example_edfs.yaml
│   ├── example_galaxy_types.yaml
│   └── example_weak_lensing.yaml
├── requirements.txt         # Python dependencies
├── README.md               # Full documentation
├── QUICKSTART.md           # Quick start guide
├── PIPELINE_OVERVIEW.md    # Technical overview
└── INSTALL.md             # This file
```

Runtime directories (created automatically):
```
├── jobs/                   # SLURM job scripts (generated)
│   ├── scripts/           # Batch scripts
│   ├── specs/             # Job specifications
│   └── logs/              # Job output logs
└── outputs/               # Results (pickle files)
```

## Troubleshooting

### TreeCorr Installation Issues

TreeCorr requires compilation. If installation fails:

**Try pre-built wheel:**
```bash
pip install treecorr --no-cache-dir
```

**Install build dependencies:**
```bash
# Ubuntu/Debian
sudo apt-get install build-essential python3-dev

# RHEL/CentOS
sudo yum install gcc gcc-c++ python3-devel

# Then retry
pip install treecorr
```

**Use conda (usually easiest):**
```bash
conda install -c conda-forge treecorr
```

### Import Errors

If you see "No module named 'src'":

```bash
# Make sure you're in the pipeline directory
cd cluster_clustering_pipeline

# Set PYTHONPATH
export PYTHONPATH=$PWD:$PYTHONPATH

# Or run with explicit path
PYTHONPATH=$PWD python scripts/run_analysis.py ...
```

### Permission Errors

```bash
# Make scripts executable
chmod +x scripts/*.py
```

### SLURM Issues

**Jobs fail immediately:**
- Check SLURM logs: `cat jobs/logs/array_*.err`
- Verify module names: `module avail python`
- Test environment activation manually

**Module not found in job:**
- Ensure `conda_env` or `modules` set correctly in config
- Test: `sbatch --wrap="python -c 'import treecorr'"`

### Memory Issues

If jobs run out of memory:

```yaml
slurm:
  mem: "32G"  # Increase as needed
```

Check actual usage:
```bash
sacct -j JOBID --format=JobID,MaxRSS,Elapsed
```

## Updating

To update the pipeline code:

```bash
cd cluster_clustering_pipeline
git pull  # if using git

# Re-check setup
python scripts/check_setup.py
```

Dependencies rarely change, but if needed:
```bash
pip install -r requirements.txt --upgrade
```

## Uninstallation

```bash
# Remove conda environment
conda env remove -n cluster_pipeline

# Or remove pip packages
pip uninstall astropy numpy treecorr pyyaml matplotlib

# Remove pipeline directory
rm -rf cluster_clustering_pipeline
```

## Next Steps

Once installation is complete:

1. Read [QUICKSTART.md](QUICKSTART.md) for a 5-minute tutorial
2. Review example configs in `configs/`
3. Create your own config: `python scripts/create_config.py`
4. Run a test: `python scripts/run_analysis.py --config config.yaml --test`
5. Generate jobs: `python scripts/generate_jobs.py --config config.yaml --array`
6. Submit: `sbatch jobs/submit_array.sh`

For detailed usage, see [README.md](README.md).

## Getting Help

- Check [QUICKSTART.md](QUICKSTART.md) for common workflows
- See [README.md](README.md) Troubleshooting section
- Review example configs in `configs/`
- Run `python scripts/check_setup.py` to diagnose issues
