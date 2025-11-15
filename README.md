# Cluster Clustering Pipeline

A modular Python pipeline for computing weighted clustering around clusters using galaxy catalogues. This pipeline supports filtering on both cluster and galaxy properties, computing clustering measurements with systematic weights, and parallel batch job submission via SLURM.

## Features

- **Flexible Filtering**: Filter clusters and galaxies by any property (redshift, richness, Sersic index, etc.)
- **Weighted Clustering**: Compute weighted cross-correlations using TreeCorr
- **Systematic Weights**: Apply spatial, depth-based, or property-based weights
- **Batch Processing**: Generate and submit SLURM job arrays for parallel processing
- **Multiple Analyses**: Support for clustering, weak lensing, and intrinsic alignment measurements
- **Configuration-Driven**: All parameters specified in YAML configuration files

## Installation

### Requirements

- Python 3.8+
- astropy
- numpy
- treecorr
- pyyaml
- matplotlib (for plotting)

### Setup

```bash
# Clone or navigate to the pipeline directory
cd cluster_clustering_pipeline

# Install dependencies (if not already installed)
pip install astropy numpy treecorr pyyaml matplotlib

# Or create a conda environment
conda create -n cluster_pipeline python=3.11
conda activate cluster_pipeline
pip install astropy numpy treecorr pyyaml matplotlib
```

## Quick Start

### 1. Create a Configuration File

```bash
python scripts/create_config.py --output my_config.yaml
```

Edit `my_config.yaml` with your catalogue paths and analysis parameters.

### 2. Test the Pipeline

Run a single job to verify everything works:

```bash
python scripts/run_analysis.py --config my_config.yaml --test
```

### 3. Generate SLURM Jobs

Create job scripts for all bin combinations:

```bash
python scripts/generate_jobs.py --config my_config.yaml --array
```

This creates:
- `jobs/scripts/` - Individual job scripts (if not using `--array`)
- `jobs/specs/` - Job specification YAML files
- `jobs/submit_array.sh` - SLURM array submission script

### 4. Submit Jobs

For job arrays:
```bash
sbatch jobs/submit_array.sh
```

For individual jobs:
```bash
bash jobs/scripts/submit_all.sh
```

### 5. Plot Results

```bash
python scripts/plot_results.py outputs/*.pkl --output clustering_plot.png
```

## Configuration File

The configuration file (YAML format) specifies all analysis parameters:

### Basic Structure

```yaml
catalogues:
  clusters: /path/to/cluster_catalogue.fits
  galaxies: /path/to/galaxy_catalogue.fits
  randoms: /path/to/random_catalogue.fits

column_mappings:
  cluster:
    ra: RIGHT_ASCENSION_CLUSTER_pzwav
    dec: DECLINATION_CLUSTER_pzwav
    redshift: Z_CLUSTER_pzwav
  galaxy:
    ra: right_ascension
    dec: declination
    redshift: phz_median

cosmology:
  H0: 70
  Om0: 0.3

analysis_parameters:
  min_sep: 0.1           # Mpc/h (or degrees if mode=2d)
  max_sep: 10.0
  nbins: 16
  min_rpar: -60          # Line-of-sight range
  max_rpar: 60
  bin_type: Log          # or Linear
  metric: Rperp          # Rperp for 3D, Euclidean for 2D
  mode: 3d               # 3d or 2d

cluster_bins:
  redshift:
    column: Z_CLUSTER_pzwav
    bins: [0.2, 0.5, 0.8, 1.2]
  richness:
    column: RICHNESS_CLUSTER_pzwav
    bins: [10, 20, 50]

galaxy_bins:
  redshift:
    column: redshift
    bins: [0.0, 0.5, 1.0, 1.5, 2.0]

output_directory: ./outputs

slurm:
  time: "02:00:00"
  mem: "16G"
  cpus: 4
  partition: "batch"
  conda_env: my_env       # Optional
```

### Binning Options

The pipeline generates jobs for **all combinations** of cluster and galaxy bins.

Example: If you specify 3 cluster redshift bins and 4 galaxy redshift bins, you get 3 × 4 = 12 jobs.

**Redshift bins:**
```yaml
cluster_bins:
  redshift:
    column: Z_CLUSTER_pzwav
    bins: [0.2, 0.5, 0.8]  # Creates 2 bins: [0.2-0.5], [0.5-0.8]
```

**Multiple properties:**
```yaml
cluster_bins:
  redshift:
    column: Z_CLUSTER_pzwav
    bins: [0.2, 0.5, 0.8]
  richness:
    column: RICHNESS_CLUSTER_pzwav
    bins: [10, 20, 50]
# Creates 2 × 2 = 4 cluster bin combinations
```

**Galaxy types (early vs late):**
```yaml
galaxy_bins:
  galaxy_type:
    column: sersic_sersic_vis_index
    bins: [0.0, 2.0, 10.0]  # Late-type: [0-2], Early-type: [2-10]
```

## Usage Examples

### Example 1: Basic Clustering

Compute clustering around clusters in different redshift bins:

```yaml
cluster_bins:
  redshift:
    column: Z_CLUSTER_pzwav
    bins: [0.2, 0.5, 0.8, 1.2]

galaxy_bins:
  redshift:
    column: redshift
    bins: [0.0, 1.0, 2.0]
```

### Example 2: Galaxy Type Dependence

Compare early-type vs late-type galaxy clustering:

```yaml
galaxy_bins:
  redshift:
    column: redshift
    bins: [0.5, 1.0, 1.5]
  galaxy_type:
    column: sersic_sersic_vis_index
    bins: [0.0, 2.0, 10.0]  # Split by Sersic index
```

### Example 3: Weak Lensing

Measure tangential shear around clusters:

```yaml
analysis_parameters:
  min_sep: 0.01          # degrees
  max_sep: 0.5
  metric: Euclidean
  sep_units: deg
  mode: 2d               # Use angular coordinates

cluster_bins:
  redshift:
    column: Z_CLUSTER_pzwav
    bins: [0.2, 0.5, 0.8]  # Lens redshifts

galaxy_bins:
  redshift:
    column: redshift
    bins: [0.8, 1.5, 2.0]  # Source redshifts (> lens z)
```

### Example 4: Systematic Weights

Apply systematic weights to correct for observational effects:

```yaml
systematic_weights:
  enabled: true
  depth_weights:
    column: limiting_magnitude
  spatial_weights:
    nside: 64
```

## Output Files

Results are saved as pickle files in the output directory:

```
outputs/
  clustering_c_Z_CLUSTER_0.20_0.50_g_redshift_0.00_1.00.pkl
  clustering_c_Z_CLUSTER_0.50_0.80_g_redshift_0.00_1.00.pkl
  ...
```

Each file contains:
- `results`: Dictionary with `r`, `xi`, `var_xi`, `npairs`
- `parameters`: Analysis parameters used
- `metadata`: Bin definitions and object counts

### Loading Results

```python
import pickle

with open('outputs/clustering_c_....pkl', 'rb') as f:
    data = pickle.load(f)

r = data['results']['r']           # Separation bins
xi = data['results']['xi']         # Correlation function
sigma = data['results']['sigma_xi'] # Uncertainties
```

## Scripts Reference

### `create_config.py`
Create a template configuration file.

```bash
python scripts/create_config.py --output my_config.yaml
```

### `run_analysis.py`
Run the full analysis locally (for testing or small jobs).

```bash
# Test with first job only
python scripts/run_analysis.py --config config.yaml --test

# Run all jobs locally
python scripts/run_analysis.py --config config.yaml

# Dry run (print what would be done)
python scripts/run_analysis.py --config config.yaml --dry-run
```

### `generate_jobs.py`
Generate SLURM job scripts.

```bash
# Generate job array (recommended for many jobs)
python scripts/generate_jobs.py --config config.yaml --array

# Generate individual job scripts
python scripts/generate_jobs.py --config config.yaml

# Specify output directory
python scripts/generate_jobs.py --config config.yaml --output-dir /path/to/jobs
```

### `plot_results.py`
Plot clustering results.

```bash
# Plot multiple results
python scripts/plot_results.py outputs/*.pkl --output plot.png

# Comparison plot
python scripts/plot_results.py file1.pkl file2.pkl --compare

# Custom labels
python scripts/plot_results.py file1.pkl file2.pkl --labels "Label 1" "Label 2"

# Linear scale
python scripts/plot_results.py outputs/*.pkl --linear
```

## Pipeline Architecture

```
cluster_clustering_pipeline/
├── src/                     # Core modules
│   ├── catalogue.py         # Load and manage catalogues
│   ├── filters.py           # Filter catalogues by properties
│   ├── clustering.py        # Compute clustering measurements
│   ├── shapes.py            # Shape-position correlations
│   ├── weights.py           # Systematic weight calculations
│   └── config.py            # Configuration management
├── scripts/                 # Command-line scripts
│   ├── create_config.py     # Create template config
│   ├── run_analysis.py      # Run full analysis
│   ├── run_single_job.py    # Run single job (called by SLURM)
│   ├── generate_jobs.py     # Generate SLURM scripts
│   └── plot_results.py      # Plot results
├── configs/                 # Example configurations
│   ├── example_edfs.yaml
│   ├── example_galaxy_types.yaml
│   └── example_weak_lensing.yaml
├── jobs/                    # Generated job scripts (created at runtime)
│   ├── scripts/
│   ├── specs/
│   └── logs/
└── outputs/                 # Results (created at runtime)
```

## Advanced Usage

### Custom Filters

You can use lambda functions for complex filtering logic:

```python
# In Python, not YAML (for now, YAML only supports ranges)
from src.filters import CatalogueFilter

custom_filters = {
    'color': lambda x: (x > 0.5) & (x < 1.5),
    'magnitude': lambda x: x < 24.0
}

filtered = CatalogueFilter.apply_filters(catalogue, custom_filters)
```

### Computing Weights Separately

```python
from src.weights import WeightCalculator, apply_systematic_weights

# Compute weights
weight_config = {
    'depth_weights': {'column': 'limiting_magnitude'},
    'spatial_weights': {'nside': 64}
}

weights = apply_systematic_weights(galaxy_catalogue, weight_config)

# Save weights for later use
calculator = WeightCalculator()
calculator.save_weights(weights, 'systematic_weights.npy',
                       metadata={'config': weight_config})
```

### Programmatic Usage

```python
from src.catalogue import CatalogueManager
from src.clustering import ClusteringAnalysis
from astropy.cosmology import FlatLambdaCDM

# Setup
cosmology = FlatLambdaCDM(H0=70, Om0=0.3)
manager = CatalogueManager(cosmology=cosmology)

# Load catalogues
clusters = manager.load_cluster_catalogue('clusters.fits')
galaxies = manager.load_galaxy_catalogue('galaxies.fits')

# Add coordinates
clusters = manager.add_cartesian_coordinates(clusters)
galaxies = manager.add_cartesian_coordinates(galaxies)

# Create TreeCorr catalogues
cluster_cat = manager.create_treecorr_catalogue(clusters, mode='3d')
galaxy_cat = manager.create_treecorr_catalogue(galaxies, mode='3d')

# Run analysis
analysis = ClusteringAnalysis(min_sep=0.1, max_sep=10, nbins=16,
                              min_rpar=-60, max_rpar=60)
results = analysis.compute_clustering(cluster_cat, galaxy_cat, random_cat)

# Save
analysis.save_results(results, 'my_results.pkl')
```

## Troubleshooting

### Issue: Jobs fail with memory errors
**Solution**: Increase `mem` in SLURM config or reduce catalogue size for testing.

### Issue: No results in output files
**Solution**: Check that filter ranges overlap with your data. Use `--test` mode to debug.

### Issue: TreeCorr warnings about patches
**Solution**: This is normal if some spatial patches are empty. The analysis will still work.

### Issue: NaN values in coordinates
**Solution**: The pipeline handles NaN redshifts gracefully by filtering them out.

## Citation

If you use this pipeline in your research, please cite:
- TreeCorr: Jarvis et al. (2004) https://github.com/rmjarvis/TreeCorr
- Your cluster detection method (e.g., PZWav, AMICO)

## Contact

For questions or issues, please contact [Your Name/Team].

## License

[Specify your license here]
