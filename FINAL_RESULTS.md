# Cluster-Galaxy Angular Clustering Analysis - Final Results

**Date**: 2025-11-15
**Analysis Type**: 2D Angular Correlation Functions with Matched Tomographic Redshift Binning
**Total Jobs**: 25 (all completed successfully)

## Analysis Configuration

### Catalogues
- **Clusters**: `/sps/euclid/Users/cmurray/RR2_alt/amico_merged_cl2025-07-17T09_36_11.fits`
  - Total: 25,843 clusters
- **Galaxies**: `/sps/euclid/Users/cmurray/rr2_data/galaxies.fits`
  - Total: 62,343,553 galaxies
- **Randoms**: `/sps/euclid/Users/cmurray/rr2_data/full_area_randoms_20250724.fits`
  - Total: 62,343,553 randoms (filtered by galaxy z-bin per job)

### Analysis Parameters
- **Mode**: 2D angular correlations
- **Metric**: Euclidean (flat-sky approximation)
- **Separation Range**: 0.1 - 60.0 arcminutes
- **Number of Bins**: 16 (logarithmic spacing)
- **Cosmology**: Flat ΛCDM with H₀=70, Ωₘ=0.3

## Binning Strategy: Matched 5×5 Tomographic Grid

**Both clusters and galaxies use identical redshift bins**:

| Bin | Redshift Range | N Clusters | N Galaxies | % Galaxies |
|-----|---------------|------------|------------|------------|
| 1   | 0.0 ≤ z < 0.3 | ~2,000     | 1,678,378  | 2.7%       |
| 2   | 0.3 ≤ z < 0.6 | ~8,000     | 13,985,587 | 22.4%      |
| 3   | 0.6 ≤ z < 0.9 | ~10,000    | 13,732,279 | 22.0%      |
| 4   | 0.9 ≤ z < 1.2 | ~4,500     | 10,364,070 | 16.6%      |
| 5   | 1.2 ≤ z < 1.5 | ~1,300     | 6,914,214  | 11.1%      |

**Total Measurements**: 25 (5 cluster bins × 5 galaxy bins)

### Correlation Matrix Structure

```
                  Galaxy Redshift Bins
              z=0.0-0.3  0.3-0.6  0.6-0.9  0.9-1.2  1.2-1.5
Cluster   0.0-0.3   [0]      [1]      [2]      [3]      [4]
Redshift  0.3-0.6   [5]      [6]      [7]      [8]      [9]
Bins      0.6-0.9  [10]     [11]     [12]     [13]     [14]
          0.9-1.2  [15]     [16]     [17]     [18]     [19]
          1.2-1.5  [20]     [21]     [22]     [23]     [24]
```

**Diagonal elements** (jobs 0, 6, 12, 18, 24): Auto-correlation within same redshift slice
**Off-diagonal elements**: Cross-redshift correlations

## Scientific Advantages of Matched Binning

1. **Auto vs Cross Correlations**:
   - Diagonal: Clusters and galaxies at same redshift (physical associations)
   - Off-diagonal: Potential contamination from photo-z errors or projection effects

2. **Photo-z Validation**:
   - Strong signal on diagonal with weak off-diagonal → good photo-z
   - Significant off-diagonal signal → photo-z scatter or catastrophic failures

3. **Redshift Evolution**:
   - Compare diagonal elements to see how clustering strength evolves with z
   - Track cluster-galaxy correlation amplitude vs redshift

4. **Systematic Tests**:
   - Off-diagonal should show expected behavior from photo-z scatter
   - Can constrain photo-z σ_z/(1+z) from off-diagonal width

## Generated Output

### Result Files
- **Location**: `outputs_rr2_angular/`
- **Format**: Pickle files (.pkl)
- **Count**: 25 files
- **Naming**: `clustering_c_Z_CLUSTER_A_B_g_redshift_C_D.pkl`

Example filenames:
- `clustering_c_Z_CLUSTER_0.00_0.30_g_redshift_0.00_0.30.pkl` (diagonal, job 0)
- `clustering_c_Z_CLUSTER_0.30_0.60_g_redshift_0.60_0.90.pkl` (off-diagonal, job 7)

Each result file contains:
```python
{
    'results': {
        'r': array([...]),           # Mean separation [arcmin]
        'xi': array([...]),          # Correlation function ω(θ)
        'var_xi': array([...]),      # Variance
        'sigma_xi': array([...]),    # Standard error
        'npairs': array([...])       # Number of pairs per bin
    },
    'parameters': {...},              # Analysis config
    'metadata': {
        'cluster_filters': {...},
        'galaxy_filters': {...},
        'n_clusters': int,
        'n_galaxies': int,
        'n_randoms': int
    }
}
```

### Diagnostic Plots

**Location**: `outputs_rr2_angular/figures/`

#### 1. Tomographic Plots (5 plots, one per cluster redshift bin)

Each plot shows the correlation function for that cluster z-bin vs all 5 galaxy z-bins:

- `tomographic_z0.00_0.30.png` - Lowest redshift clusters
- `tomographic_z0.30_0.60.png`
- `tomographic_z0.60_0.90.png` - Peak cluster redshift
- `tomographic_z0.90_1.20.png`
- `tomographic_z1.20_1.50.png` - Highest redshift clusters

**Features**:
- Top panel: Linear scale showing full correlation function
- Bottom panel: Log-log scale highlighting power-law behavior
- Color-coded by galaxy redshift bin (5 different colors)
- Error bars from TreeCorr jackknife variance
- Diagonal element typically strongest (same redshift)

#### 2. Summary Grid (`clustering_summary_grid.png`)

**Comprehensive 5-panel figure** showing all cluster redshift bins:
- Each panel = one cluster z-bin
- Each panel shows 5 galaxy z-bins overlaid
- Legend indicates galaxy redshift ranges
- Compact overview of entire 5×5 matrix

**How to Read**:
- Look for strongest signal on diagonal (e.g., panel for z_cl=0.6-0.9 should have strongest curve for z_gal=0.6-0.9)
- Off-diagonal strength indicates photo-z scatter
- Evolution across panels shows redshift dependence

#### 3. Galaxy N(z) Plot (`galaxy_nz_all_bins.png`)

Shows redshift distribution for all 5 galaxy tomographic bins:
- Normalized histograms with different colors per bin
- Statistics printed for each bin (mean, median, std)
- Verifies bin definitions and galaxy counts

## Performance Metrics

### Computational Efficiency
- **Jobs**: 25 (5×5 grid)
- **Time per job**: ~4-5 minutes
- **Total runtime**: ~20-30 minutes (25 jobs in parallel)
- **Memory usage**: ~50-60GB per job
- **Resources**: 8 CPUs, 64GB RAM per job

### Optimizations Applied
1. **2D Angular Mode**: No Cartesian coordinate calculation
2. **Column Loading**: Only RA, Dec, redshift (3 columns vs 100+)
3. **Random Filtering**: Randoms filtered by galaxy z-bin (~80% memory reduction)
4. **Euclidean Metric**: Flat-sky approximation (faster than Arc)

**Speed Improvement**: 25 jobs in ~25 mins vs previous 90 jobs in ~2 hours

## Usage Examples

### Load and Plot a Single Result

```python
import pickle
import matplotlib.pyplot as plt
import numpy as np

# Load diagonal element (same redshift)
with open('outputs_rr2_angular/clustering_c_Z_CLUSTER_0.60_0.90_g_redshift_0.60_0.90.pkl', 'rb') as f:
    data = pickle.load(f)

r = data['results']['r']              # Separation [arcmin]
xi = data['results']['xi']            # Correlation function
sigma_xi = data['results']['sigma_xi'] # Error

# Plot
plt.errorbar(r, xi, yerr=sigma_xi, marker='o', capsize=3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('θ [arcmin]')
plt.ylabel('ω(θ)')
plt.title('Cluster-Galaxy Clustering (z=0.6-0.9)')
plt.grid(True, alpha=0.3)
plt.show()
```

### Compare Diagonal vs Off-Diagonal

```python
import pickle
import matplotlib.pyplot as plt

# Load diagonal (z_cl = z_gal = 0.6-0.9)
with open('outputs_rr2_angular/clustering_c_Z_CLUSTER_0.60_0.90_g_redshift_0.60_0.90.pkl', 'rb') as f:
    diag = pickle.load(f)

# Load off-diagonal (z_cl = 0.6-0.9, z_gal = 0.3-0.6)
with open('outputs_rr2_angular/clustering_c_Z_CLUSTER_0.60_0.90_g_redshift_0.30_0.60.pkl', 'rb') as f:
    off_diag = pickle.load(f)

# Plot comparison
plt.errorbar(diag['results']['r'], diag['results']['xi'],
             yerr=diag['results']['sigma_xi'],
             marker='o', label='Same z (0.6-0.9)', capsize=3)
plt.errorbar(off_diag['results']['r'], off_diag['results']['xi'],
             yerr=off_diag['results']['sigma_xi'],
             marker='s', label='Cl: 0.6-0.9, Gal: 0.3-0.6', capsize=3)
plt.xscale('log')
plt.xlabel('θ [arcmin]')
plt.ylabel('ω(θ)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.title('Diagonal vs Off-Diagonal Clustering')
plt.show()
```

### Regenerate All Plots

```bash
cd /sps/euclid/Users/cmurray/euclid_cluster_clustering/cluster_clustering_pipeline

python scripts/plot_results.py \
    --results-dir outputs_rr2_angular \
    --output-dir outputs_rr2_angular/figures \
    --plot-type all
```

## Key Science Questions

With this matched tomographic analysis, you can now investigate:

1. **Redshift Evolution of Clustering**:
   - How does cluster-galaxy correlation strength change with z?
   - Compare diagonal elements: jobs 0, 6, 12, 18, 24

2. **Photo-z Quality**:
   - How strong is off-diagonal signal compared to diagonal?
   - Estimate photo-z scatter from off-diagonal "wings"

3. **Physical Associations**:
   - Are clusters preferentially surrounded by galaxies at same z?
   - Expected from collapse of overdense regions

4. **Projection Effects**:
   - How much contamination from foreground/background galaxies?
   - Off-diagonal signal at large scales = projections

5. **Cluster Bias**:
   - Fit clustering amplitude to measure cluster bias b(z)
   - Compare to theoretical predictions

## Next Steps

### Immediate Follow-ups:
1. Create 5×5 heatmap of clustering amplitude
2. Fit power-law models: ω(θ) = A θ^(-γ)
3. Extract correlation length r₀ and slope γ per bin
4. Plot diagonal vs off-diagonal amplitudes
5. Compare to theoretical predictions

### Advanced Analysis:
1. Halo model fitting to constrain M_halo
2. Cross-correlate with weak lensing
3. Split by cluster richness (add richness dimension)
4. Systematic tests (depth, seeing, photo-z quality)

## Files and Directory Structure

```
cluster_clustering_pipeline/
├── configs/
│   └── rr2_angular_tomographic.yaml  # 5×5 matched bins config
├── outputs_rr2_angular/
│   ├── clustering_c_*.pkl            # 25 result files
│   └── figures/
│       ├── galaxy_nz_all_bins.png
│       ├── clustering_summary_grid.png
│       └── tomographic/
│           ├── tomographic_z0.00_0.30.png
│           ├── tomographic_z0.30_0.60.png
│           ├── tomographic_z0.60_0.90.png
│           ├── tomographic_z0.90_1.20.png
│           └── tomographic_z1.20_1.50.png
├── scripts/
│   ├── plot_results.py               # Plotting script
│   ├── plot_nz.py                    # N(z) diagnostic
│   └── run_single_job.py             # Job execution
└── jobs/
    ├── submit_array.sh               # SLURM submission script
    ├── specs/                        # 25 job specifications
    └── logs/                         # Job output logs
```

## Summary

✅ **25 measurements completed** in 5×5 matched tomographic grid
✅ **All plots generated**: 5 tomographic + 1 summary grid + 1 n(z)
✅ **Memory optimized**: 64GB sufficient with column loading
✅ **Fast execution**: ~25 minutes total runtime

This matched binning approach provides a clean framework for studying:
- Redshift evolution of cluster-galaxy clustering
- Photo-z quality and contamination
- Physical associations vs projection effects
- Systematic validation of clustering measurements

The results are ready for scientific analysis and publication!
