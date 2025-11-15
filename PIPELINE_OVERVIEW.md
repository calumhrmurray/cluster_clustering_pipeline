# Pipeline Overview

## What This Pipeline Does

This pipeline automates the computation of **weighted clustering measurements** around galaxy clusters. It:

1. **Loads** cluster and galaxy catalogues from FITS files
2. **Filters** objects based on properties (redshift, richness, galaxy type, etc.)
3. **Computes** cross-correlation functions using TreeCorr
4. **Applies** systematic weights to correct for observational effects
5. **Parallelizes** analysis across all bin combinations via SLURM
6. **Saves** results for later analysis and plotting

## Key Concepts

### Clustering Measurement

The pipeline computes the cluster-galaxy cross-correlation function using the **Landy-Szalay estimator**:

```
ξ(r) = [DD - DR - RD + RR] / RR
```

Where:
- **DD** = Data-Data pairs (clusters × galaxies)
- **DR** = Data-Random pairs (clusters × randoms)
- **RD** = Random-Data pairs (randoms × galaxies)
- **RR** = Random-Random pairs

### Binning Strategy

The pipeline creates **independent jobs** for each combination of cluster and galaxy bins:

```
Cluster bins:  [z: 0.2-0.5, z: 0.5-0.8] × [rich: 10-20, rich: 20-50]
Galaxy bins:   [z: 0-1, z: 1-2]

Total jobs = 2 × 2 × 2 = 8
```

Each job:
1. Filters catalogues to the specified bin
2. Computes clustering for that specific combination
3. Saves results independently

This allows:
- **Parallel processing** on HPC clusters
- **Easy restart** of failed jobs
- **Flexible post-processing** (combine bins as needed)

### Coordinate Systems

**3D Mode (Cartesian coordinates)**:
- Converts RA, Dec, z → x, y, z in Mpc/h
- Uses `metric='Rperp'` for projected clustering
- Separates line-of-sight (r_parallel) from transverse (r_perp) distances
- Best for: Clustering analysis with redshift information

**2D Mode (Angular coordinates)**:
- Uses RA, Dec directly in degrees
- Uses `metric='Euclidean'` with angular separations
- Projection along line of sight
- Best for: Weak lensing, where projection effects matter

## Data Flow

```
Input FITS Catalogues
        ↓
    Load & Standardize
        ↓
    Apply Filters (by bin)
        ↓
    Add Coordinates (Cartesian or Angular)
        ↓
    Apply Systematic Weights (optional)
        ↓
    Create TreeCorr Catalogues
        ↓
    Compute Correlations (DD, DR, RD, RR)
        ↓
    Calculate ξ(r) and uncertainties
        ↓
    Save Results (pickle files)
        ↓
    Plot & Analyze
```

## Module Responsibilities

### `catalogue.py` - Data Management
- Load FITS catalogues
- Standardize column names
- Add Cartesian/angular coordinates
- Convert ellipticities to shear convention
- Create TreeCorr catalogue objects

### `filters.py` - Selection Functions
- Apply range filters: `property: (min, max)`
- Apply categorical filters: `property: [value1, value2]`
- Select galaxy types (early/late)
- Create lens-source splits
- Generate human-readable bin labels

### `clustering.py` - Correlation Analysis
- Compute NN (number-number) correlations
- Implement Landy-Szalay estimator
- Handle weighted correlations
- Save/load results with metadata

### `shapes.py` - Lensing Measurements
- Compute NG (number-shear) correlations
- Measure tangential shear profiles
- Intrinsic alignment analysis
- Support for cross/plus components

### `weights.py` - Systematic Corrections
- Depth-based weights
- Spatial uniformity weights
- Property distribution matching
- Weight combination and normalization

### `config.py` - Configuration Management
- Parse YAML configuration files
- Generate job specifications
- Validate configuration
- Create bin combinations

## File Naming Convention

Results are saved with descriptive filenames encoding the bin parameters:

```
clustering_c_{cluster_label}_g_{galaxy_label}.pkl
```

Example:
```
clustering_c_Z_CLUSTER_0.20_0.50_RICHNESS_10.00_20.00_g_redshift_0.50_1.00.pkl
```

This means:
- Clusters: 0.2 < z < 0.5, 10 < richness < 20
- Galaxies: 0.5 < z < 1.0

## SLURM Integration

### Job Array Mode (Recommended)

Single submission script runs all jobs:
```bash
sbatch jobs/submit_array.sh
```

Each array task:
- Reads its job spec from `jobs/specs/job_N_spec.yaml`
- Runs `run_single_job.py` with that spec
- Independent of other jobs (no dependencies)

**Advantages**:
- One submission for all jobs
- Easy to manage
- Automatic load balancing
- Can resubmit failed jobs individually

### Individual Job Mode

Generate separate scripts for each job:
```bash
python scripts/generate_jobs.py --config config.yaml
bash jobs/scripts/submit_all.sh
```

**Use when**:
- Jobs have different resource requirements
- Need custom dependencies
- Want fine-grained control

## Extending the Pipeline

### Adding New Filter Types

In `filters.py`, add a method to `CatalogueFilter`:

```python
@staticmethod
def select_color_cut(catalogue, color_col='color', threshold=0.5):
    """Select galaxies by color."""
    filters = {color_col: lambda x: x > threshold}
    return CatalogueFilter.apply_filters(catalogue, filters)
```

### Adding New Weight Types

In `weights.py`, add a method to `WeightCalculator`:

```python
def compute_color_weights(self, catalogue, color_col='color',
                         reference_distribution=None):
    """Weight by color distribution."""
    # Your implementation
    return weights
```

### Custom Analysis Types

Create a new module (e.g., `src/my_analysis.py`) following the pattern:

```python
class MyAnalysis:
    def __init__(self, **params):
        self.params = params

    def compute(self, cat1, cat2):
        # Your analysis
        return results

    def save_results(self, results, path):
        # Save implementation
        pass
```

## Performance Considerations

### Memory Usage

Each job loads:
- Full cluster catalogue
- Full galaxy catalogue
- Full random catalogue

Then filters to the relevant bin. For very large catalogues:

**Option 1**: Pre-filter catalogues
```bash
# Create separate files per redshift slice
python scripts/prefilter_catalogues.py
```

**Option 2**: Increase memory allocation
```yaml
slurm:
  mem: "32G"  # or higher
```

### Computation Time

Typical scaling:
- Small job (1k clusters, 100k galaxies): ~5-10 minutes
- Medium job (10k clusters, 1M galaxies): ~30-60 minutes
- Large job (100k clusters, 10M galaxies): ~2-4 hours

Factors:
- Number of bins (nbins)
- Separation range
- Number of patches (for jackknife)
- 3D vs 2D mode

### Parallelization

Within each job, TreeCorr uses multi-threading:
```yaml
slurm:
  cpus: 8  # TreeCorr will use 8 threads
```

Across jobs, SLURM handles parallelization:
```bash
# 100 jobs, run 20 at a time
#SBATCH --array=0-99%20
```

## Output Analysis

### Combining Results

To stack results from multiple bins:

```python
import pickle
import numpy as np
from pathlib import Path

results = []
for file in Path('outputs').glob('clustering_*.pkl'):
    with open(file, 'rb') as f:
        results.append(pickle.load(f))

# Extract and average
all_xi = [r['results']['xi'] for r in results]
mean_xi = np.mean(all_xi, axis=0)
```

### Statistical Analysis

Results include jackknife covariance from TreeCorr patches:

```python
data = pickle.load(open('result.pkl', 'rb'))
r = data['results']['r']
xi = data['results']['xi']
cov = data['results']['var_xi']  # Variance from jackknife

# Significance
significance = xi / np.sqrt(cov)
```

### Comparison Plots

Use the provided plotting script:

```bash
# Compare different redshift bins
python scripts/plot_results.py \
  outputs/*_z_0.2_0.5*.pkl \
  outputs/*_z_0.5_0.8*.pkl \
  --labels "z=0.2-0.5" "z=0.5-0.8" \
  --compare
```

## Future Development

Planned features:
- [ ] HEALPix-based spatial weights (in `weights.py`)
- [ ] Automatic covariance matrix computation
- [ ] Support for FITS table output (in addition to pickle)
- [ ] Interactive plotting dashboard
- [ ] Profile stacking across redshift bins
- [ ] Support for additional correlation types (GG, NG with IA)

## References

- **TreeCorr**: Jarvis et al. (2004), https://github.com/rmjarvis/TreeCorr
- **Landy-Szalay Estimator**: Landy & Szalay (1993), ApJ 412, 64
- **Cluster Detection**: Euclid Collaboration papers on AMICO, PZWav methods
