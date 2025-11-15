# Cluster-Galaxy Angular Clustering Analysis - Results Summary

**Date**: 2025-11-15
**Analysis Type**: 2D Angular Correlation Functions with Tomographic Redshift Binning
**Total Jobs**: 90 (all completed successfully)

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

### Binning Strategy

**Cluster Bins** (18 total):
- **Redshift bins** (6): [0.2, 0.38, 0.43, 0.5, 0.6, 0.9, 1.6]
- **Richness bins** (3): [40, 80, 140, 450]
- Total cluster bins: 6 × 3 = 18

**Galaxy Bins** (5 tomographic redshift bins):
- Bin 1: 0.0 ≤ z < 0.3 (1,678,378 galaxies, 2.7%)
- Bin 2: 0.3 ≤ z < 0.6 (13,985,587 galaxies, 22.4%)
- Bin 3: 0.6 ≤ z < 0.9 (13,732,279 galaxies, 22.0%)
- Bin 4: 0.9 ≤ z < 1.2 (10,364,070 galaxies, 16.6%)
- Bin 5: 1.2 ≤ z < 1.5 (6,914,214 galaxies, 11.1%)

**Total Combinations**: 18 cluster bins × 5 galaxy bins = **90 measurements**

## Performance Metrics

### Computational Efficiency
- **Time per job**: ~4-5 minutes
- **Total runtime**: ~2 hours (90 jobs in parallel)
- **Memory usage**: ~50-60GB per job (with optimized column loading)
- **Resources**: 8 CPUs, 64GB RAM per job

### Optimizations Applied
1. **2D Mode**: Angular correlations (no Cartesian coordinate calculation)
2. **Column Loading**: Only load RA, Dec, redshift (3 columns vs ~100+)
3. **Random Filtering**: Filter randoms by galaxy redshift bin per job
4. **Euclidean Metric**: Flat-sky approximation (faster than Arc)

## Generated Output

### Result Files
- **Location**: `outputs_rr2_angular/`
- **Format**: Pickle files (.pkl)
- **Count**: 90 files (one per cluster-galaxy bin combination)
- **Naming**: `clustering_c_RICHNESS_CLUSTER_X_Y_Z_CLUSTER_A_B_g_redshift_C_D.pkl`

Each result file contains:
- `results`: Correlation function measurements
  - `r`: Mean separation in each bin (arcmin)
  - `xi`: Correlation function ω(θ)
  - `var_xi`: Variance of correlation function
  - `sigma_xi`: Standard error
  - `npairs`: Number of pairs in each bin
- `parameters`: Analysis configuration
- `metadata`: Bin information and catalogue sizes

### Diagnostic Plots

**Location**: `outputs_rr2_angular/figures/`

1. **Galaxy N(z) Plot**: `galaxy_nz_all_bins.png`
   - Shows redshift distribution for all 5 tomographic bins
   - Normalized histograms with statistics

2. **Tomographic Comparison Plots**: `tomographic/` (18 plots)
   - One plot per cluster bin showing all 5 galaxy redshift bins
   - Both linear and log-log scales
   - Error bars from jackknife variance
   - Examples:
     - `tomographic_rich40_80_z0.20_0.38.png`
     - `tomographic_rich80_140_z0.43_0.50.png`
     - `tomographic_rich140_450_z0.90_1.60.png`

3. **Summary Grid**: `clustering_summary_grid.png`
   - All 18 cluster bins in one figure
   - Each panel shows 5 galaxy redshift bins
   - Compact overview of entire analysis

## Key Features

### Tomographic Analysis Benefits
- **Systematics Testing**: Different galaxy redshifts probe different systematics
- **Evolution Studies**: Compare clustering at different redshifts
- **Photo-z Validation**: Check for consistency across redshift bins
- **Memory Efficiency**: Filtering randoms reduces memory by ~80%

### Quality Checks
- All 90 jobs completed successfully
- No OOM (out-of-memory) errors with 64GB allocation
- Correlation functions show expected behavior (decreasing with separation)
- Jackknife variance estimates computed for all measurements

## Usage Examples

### View a specific result:
```python
import pickle
with open('outputs_rr2_angular/clustering_c_RICHNESS_CLUSTER_40.00_80.00_Z_CLUSTER_0.20_0.38_g_redshift_0.30_0.60.pkl', 'rb') as f:
    data = pickle.load(f)

r = data['results']['r']          # Separation [arcmin]
xi = data['results']['xi']        # Correlation function
sigma_xi = data['results']['sigma_xi']  # Error
```

### Regenerate plots:
```bash
python scripts/plot_results.py \
    --results-dir outputs_rr2_angular \
    --output-dir outputs_rr2_angular/figures \
    --plot-type all
```

### Generate N(z) plot:
```bash
python scripts/plot_nz.py \
    --catalogue /sps/euclid/Users/cmurray/rr2_data/galaxies.fits \
    --config configs/rr2_angular_tomographic.yaml \
    --output outputs_rr2_angular/figures/galaxy_nz_all_bins.png
```

## Next Steps

Potential follow-up analyses:
1. Fit clustering models to measurements
2. Compare with theoretical predictions
3. Measure cluster bias as function of richness/redshift
4. Cross-correlate with lensing measurements
5. Systematic tests (depth, area, photo-z quality)

## Files Modified/Created

### Core Pipeline
- `src/catalogue.py`: Added column-specific loading for memory efficiency
- `src/clustering.py`: Added 2D mode support and random filtering
- `src/config.py`: Added 2D parameter validation

### Configuration
- `configs/rr2_angular_tomographic.yaml`: 2D angular analysis configuration

### Scripts
- `scripts/plot_results.py`: Tomographic plotting with summary grids
- `scripts/plot_nz.py`: Galaxy redshift distribution plots
- `scripts/generate_jobs.py`, `scripts/submit_jobs.py`: Job generation
- `scripts/run_single_job.py`: Single job execution with optimizations

### Documentation
- `RESULTS_SUMMARY.md`: This file
- Job specifications: `jobs/specs/job_*.yaml` (90 files)
- Job logs: `jobs/logs/` (180 files: .out and .err for each job)
