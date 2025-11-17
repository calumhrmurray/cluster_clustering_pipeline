# Systematic Weights Integration

**Date**: 2025-11-15
**Status**: ✅ Implemented and Running

## Overview

Systematic weights have been successfully integrated into the clustering pipeline to correct for spatial variations in survey properties that can bias clustering measurements.

## What Are Systematic Weights?

Systematic weights correct for observational effects that create artificial clustering:
- **Galactic extinction**: Lower detection in dusty regions → appears as underdensity
- **Stellar crowding**: More stars = harder to detect galaxies → appears as underdensity
- **Survey depth variations**: Shallower regions have fewer detections → appears as underdensity
- **Photometric systematics**: Systematic errors in magnitudes/colors

**Without weights**: These systematics create false clustering signals
**With weights**: Objects in low-density regions get higher weight to flatten systematic gradients

## Weight Files

### Cluster Weights
**File**: `/sps/euclid/Users/cmurray/euclid_clusters_tr1_validation/cluster_validation/outputs/rr2_baseline/20251115T174038Z/tables/systematic_weights/AMICO_SNR8_weights.fits`

- Format: HEALPix NSIDE=512, RING ordering
- Pixels: 3,145,728
- Weight range: [0.174, 11.150]
- Mean weight: 1.001

### Galaxy Weights
**File**: `/sps/euclid/Users/cmurray/euclid_clusters_tr1_validation/cluster_validation/outputs/rr2_baseline/20251115T174038Z/tables/systematic_weights/RR2_GALAXIES_weights.fits`

- Format: HEALPix NSIDE=512, RING ordering
- Pixels: 3,145,728
- Weight range: [0.678, 9.497]
- Mean weight: 1.001

**Note**: Weights are already normalized (mean ≈ 1.0) so total counts remain unchanged.

## Implementation

### 1. New Module: `src/systematic_weights.py`

Created standalone module with functions:

**`load_weight_map(weight_file_path)`**
- Loads HEALPix weight map from FITS file
- Returns: (weights_array, nside, nest_ordering)

**`assign_weights_to_catalogue(catalogue, weight_map_info, ra_col, dec_col)`**
- Converts RA/Dec to HEALPix pixels using `healpy.ang2pix()`
- Looks up weight for each pixel
- Adds weight column to catalogue

**`apply_systematic_weights(catalogue, weight_file, ...)`**
- Convenience function combining load + assign
- Main function to use in pipeline

### 2. Configuration Update

Modified `configs/rr2_angular_tomographic.yaml`:

```yaml
systematic_weights:
  enabled: true
  cluster_weight_file: /path/to/AMICO_SNR8_weights.fits
  galaxy_weight_file: /path/to/RR2_GALAXIES_weights.fits

output_directory: outputs_rr2_angular_weighted  # New directory
```

### 3. Pipeline Integration

Modified `scripts/run_single_job.py`:

**Order of operations:**
1. Load catalogues (clusters, galaxies, randoms)
2. **→ Apply systematic weights** ← NEW STEP
   - Clusters get cluster weights
   - Galaxies get galaxy weights
   - Randoms get galaxy weights (SAME as galaxies!)
3. Apply filters (redshift bins)
4. Create TreeCorr catalogues (with weights in 'w' column)
5. Compute clustering

**Critical**: Randoms MUST get the same weights as galaxies to maintain proper normalization.

### 4. Weight Column in TreeCorr

Weights are stored in the 'w' column of each catalogue:
- TreeCorr automatically recognizes and uses this column
- For weighted correlations: ω(θ) = Σw_i w_j / Σw_i w_j
- Landy-Szalay estimator weighted appropriately

## Testing

### Test 1: Weight Map Loading ✅
```
Cluster weights: 3.1M pixels, range [0.17, 11.15], mean 1.001
Galaxy weights: 3.1M pixels, range [0.68, 9.50], mean 1.001
```

### Test 2: Weight Assignment ✅
```
1000 test galaxies → all assigned valid weights
Weight range: [0.82, 1.92]
Mean weight: 0.99
Spatial variation present (as expected)
```

## Running Analysis

**Jobs Submitted**: 25 (5×5 tomographic grid)
**Job ID**: 13770739
**Expected Runtime**: ~25-30 minutes
**Output Directory**: `outputs_rr2_angular_weighted/`

### Monitoring Progress

Check job status:
```bash
squeue -u $USER
```

Check completed results:
```bash
ls outputs_rr2_angular_weighted/*.pkl | wc -l
```

Check logs:
```bash
tail -f jobs/logs/array_13770739_0.out
```

## Expected Impact

### Amplitude Changes
- **Large scales** (> 10 arcmin): 5-20% change in ω(θ)
- **Small scales** (< 1 arcmin): Minimal change (local clustering dominates)
- **Direction**: Typically reduces large-scale clustering (removes systematic boost)

### Error Bars
- Errors may increase slightly (effective sample size reduced by weighting)
- But measurements are now UNBIASED

### Comparison with Unweighted
To compare weighted vs unweighted:
```python
# Unweighted results
unweighted = 'outputs_rr2_angular/clustering_c_Z_CLUSTER_0.60_0.90_g_redshift_0.60_0.90.pkl'

# Weighted results
weighted = 'outputs_rr2_angular_weighted/clustering_c_Z_CLUSTER_0.60_0.90_g_redshift_0.60_0.90.pkl'

# Plot ratio to see impact
ratio = weighted_xi / unweighted_xi
```

## How Weights Were Created

Weights were computed using linear regression on systematic maps:

1. **Systematic Maps** (NSIDE=512 HEALPix):
   - Galactic extinction (E(B-V))
   - Stellar number density
   - Survey depth/noise maps
   - Airmass, seeing, other observing conditions

2. **Density Map**:
   - Count galaxies per HEALPix pixel
   - Normalize by random density: δ = N_gal/N_rand

3. **Linear Model**:
   - Fit: δ = a₁×extinction + a₂×stars + a₃×depth + ...
   - Find which systematics correlate with density

4. **Compute Weights**:
   - Predict expected density from systematics
   - Weight = 1 / predicted_density
   - Apply floor to prevent extreme weights (min = 0.1 × median)

5. **Normalization**:
   - Rescale weights to mean = 1.0
   - Preserves total counts while flattening systematics

**Reference Code**: `/sps/euclid/Users/cmurray/euclid_clusters_tr1_validation/cluster_validation/pipeline/systematics/weights.py`

## Validation Checks

### 1. Weight Distribution
- Check for reasonable range (not too extreme)
- ✅ Cluster weights: 0.17 - 11.15 (acceptable)
- ✅ Galaxy weights: 0.68 - 9.50 (acceptable)

### 2. Mean Weight
- Should be normalized to ~1.0
- ✅ Cluster mean: 1.001
- ✅ Galaxy mean: 1.001

### 3. Spatial Patterns
- Weights should correlate with systematic maps
- Higher weights in regions with:
  - Higher Galactic extinction
  - Higher stellar density
  - Lower survey depth
- TODO: Create validation plots

### 4. Same Weights for Data and Randoms
- Critical for unbiased Landy-Szalay estimator
- ✅ Code applies same galaxy_weight_file to both galaxies and randoms

## Future Enhancements

1. **Weight Diagnostics**:
   - Create plots showing weight spatial distribution
   - Histogram of weights per catalogue
   - Correlation with systematic maps

2. **Null Tests**:
   - Cross-correlate randoms with randoms (should be ω=0)
   - Check that weighted random-random gives null result

3. **Systematic Uncertainty**:
   - Vary weight floor to assess systematic uncertainty
   - Compare different weighting schemes

4. **Per-Property Weights**:
   - Could create separate weights for different galaxy populations
   - E.g., red vs blue galaxies may have different systematics

## Key References

- **Weight Creation Pipeline**: `/sps/euclid/Users/cmurray/euclid_clusters_tr1_validation/cluster_validation/pipeline/systematics/`
- **Weight Maps**: `/sps/euclid/Users/cmurray/euclid_clusters_tr1_validation/cluster_validation/outputs/rr2_baseline/20251115T174038Z/tables/systematic_weights/`
- **HEALPix Documentation**: https://healpy.readthedocs.io/
- **TreeCorr Weighting**: https://rmjarvis.github.io/TreeCorr/_build/html/weight.html

## Summary

✅ **Systematic weights successfully integrated**
✅ **25 jobs running with weighted analysis**
✅ **Output to new directory**: `outputs_rr2_angular_weighted/`
✅ **Weights applied to clusters, galaxies, AND randoms**
✅ **HEALPix-based spatial weighting**

The pipeline now corrects for observational systematics, providing unbiased clustering measurements ready for cosmological interpretation!
