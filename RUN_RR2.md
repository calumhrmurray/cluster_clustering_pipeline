# Running RR2 Cluster-Galaxy Clustering Analysis

This guide provides exact commands to run the clustering pipeline with your RR2 data.

## Quick Start (3 Commands)

```bash
cd /sps/euclid/Users/cmurray/euclid_cluster_clustering/cluster_clustering_pipeline
python scripts/run_analysis.py --config configs/rr2_clustering.yaml --test
python scripts/submit_jobs.py --config configs/rr2_clustering.yaml
```

That's it! Monitor with `squeue -u $USER`

## Data Locations

- **Clusters**: `/sps/euclid/Users/cmurray/RR2_alt/amico_merged_cl2025-07-17T09_36_11.fits`
- **Galaxies**: `/sps/euclid/Users/cmurray/rr2_data/galaxies.fits`
- **Randoms**: `/sps/euclid/Users/cmurray/rr2_data/full_area_randoms_20250724.fits`

## Configuration

Configuration file is ready at: `configs/rr2_clustering.yaml`

This creates jobs for:
- **Cluster redshift bins**: [0.2-0.38], [0.38-0.43], [0.43-0.5], [0.5-0.6], [0.6-0.9], [0.9-1.6] (6 bins)
- **Cluster richness bins**: [40-80], [80-140], [140-450] (3 bins)
- **Galaxy types**: Early-type (Sersic>2), Late-type (Sersic<2) (2 bins)

**Total jobs**: 6 × 3 × 2 = **36 jobs**

## Step-by-Step Commands

### 1. Navigate to Pipeline Directory

```bash
cd /sps/euclid/Users/cmurray/euclid_cluster_clustering/cluster_clustering_pipeline
```

### 2. Check Setup (Optional)

Verify all dependencies are installed:

```bash
python scripts/check_setup.py
```

If treecorr is missing:
```bash
pip install treecorr
```

### 3. Test with a Single Job

Run a quick test to make sure everything works:

```bash
python scripts/run_analysis.py --config configs/rr2_clustering.yaml --test
```

This will:
- Load your catalogues
- Run the first bin combination
- Save results to `outputs_rr2/`

**Check the output**: Look for a `.pkl` file in `outputs_rr2/`

### 4. Generate and Submit SLURM Jobs

Generate and submit all jobs in one command:

```bash
python scripts/submit_jobs.py --config configs/rr2_clustering.yaml
```

This will:
1. Generate 36 job specification files (one per bin combination)
2. Create SLURM array submission script
3. Submit the job array to SLURM
4. Show you the job ID and monitoring commands

**Dry run** (generate but don't submit):
```bash
python scripts/submit_jobs.py --config configs/rr2_clustering.yaml --dry-run
```

**Review what will run**:
```bash
ls jobs/specs/ | wc -l  # Should show 36
cat jobs/specs/job_0_spec.yaml  # Check first job
```

**Monitor jobs**:
```bash
# Check job status
squeue -u $USER

# Check specific job array
squeue -u $USER | grep cluster_clust

# Count running/pending jobs
squeue -u $USER -t RUNNING | wc -l
squeue -u $USER -t PENDING | wc -l
```

**Check logs** (while running or after completion):
```bash
# View first job output
tail -f jobs/logs/array_*_0.out

# Check for errors
grep -i error jobs/logs/*.err

# List completed jobs
ls jobs/logs/*.out | wc -l
```

### 5. Check Results

Once jobs complete, check the output:

```bash
# Count output files (should be 36)
ls outputs_rr2/*.pkl | wc -l

# List results
ls -lh outputs_rr2/

# Check file sizes (should be a few hundred KB each)
du -sh outputs_rr2/
```

### 6. Plot Results

Plot clustering results:

```bash
# Plot all results
python scripts/plot_results.py outputs_rr2/*.pkl --output rr2_clustering_all.png

# Compare different redshift bins (example)
python scripts/plot_results.py \
  outputs_rr2/clustering_c_Z_CLUSTER_0.20_0.38*.pkl \
  outputs_rr2/clustering_c_Z_CLUSTER_0.38_0.43*.pkl \
  outputs_rr2/clustering_c_Z_CLUSTER_0.43_0.50*.pkl \
  --compare --output rr2_redshift_comparison.png

# Compare galaxy types
python scripts/plot_results.py \
  outputs_rr2/*sersic_*_0.00_2.00.pkl \
  outputs_rr2/*sersic_*_2.00_10.00.pkl \
  --labels "Late-type" "Early-type" \
  --compare --output rr2_galaxy_types.png

# Compare richness bins
python scripts/plot_results.py \
  outputs_rr2/*RICHNESS_40*80*.pkl \
  outputs_rr2/*RICHNESS_80*140*.pkl \
  outputs_rr2/*RICHNESS_140*450*.pkl \
  --labels "λ=40-80" "λ=80-140" "λ=140-450" \
  --compare --output rr2_richness_comparison.png
```

## Understanding the Job Array

The job array works as follows:
- SLURM_ARRAY_TASK_ID ranges from 0 to 53
- Each task runs independently
- Each task reads its specification from `jobs/specs/job_N_spec.yaml`
- Results are saved to separate files

**Example job spec** (`jobs/specs/job_0_spec.yaml`):
```yaml
job_id: 0
cluster_filters:
  Z_CLUSTER: [0.2, 0.38]
  RICHNESS_CLUSTER: [40, 80]
galaxy_filters:
  sersic_sersic_vis_index: [0.0, 2.0]
```

## Troubleshooting

### Jobs fail immediately

Check SLURM logs:
```bash
cat jobs/logs/array_*_0.err
```

Common issues:
- **Module not found**: Add module loads to config SLURM section
- **Memory error**: Increase `mem` in config
- **File not found**: Check catalogue paths in config

### Some jobs succeed, others fail

Check which jobs failed:
```bash
# Check all error logs
for f in jobs/logs/array_*.err; do
  if [ -s "$f" ]; then
    echo "$f has errors:"
    tail -5 "$f"
  fi
done
```

Resubmit specific failed jobs:
```bash
# If job 5 failed, resubmit just that one:
sbatch --array=5 jobs/submit_array.sh
```

### Not enough results

Check how many jobs completed successfully:
```bash
# Count successes in logs
grep -l "JOB COMPLETED SUCCESSFULLY" jobs/logs/*.out | wc -l

# Compare with output files
ls outputs_rr2/*.pkl | wc -l
```

### Out of memory

Edit config to increase memory:
```yaml
slurm:
  mem: "64G"  # Increase from 32G
```

Regenerate jobs:
```bash
python scripts/generate_jobs.py --config configs/rr2_clustering.yaml --array
sbatch jobs/submit_array.sh
```

## Expected Runtime

Based on your catalogue sizes:
- Small bins (few clusters): ~10-20 minutes per job
- Large bins (many clusters): ~1-2 hours per job
- Total wall time (with 54 jobs in parallel): ~2-4 hours

## Modifying Bins

To change bin definitions, edit `configs/rr2_clustering.yaml`:

```yaml
# Example: Coarser redshift bins
cluster_bins:
  redshift:
    column: Z_CLUSTER
    bins: [0.2, 0.5, 0.9, 1.6]  # 3 bins instead of 6

# Example: No richness binning (use all clusters)
cluster_bins:
  redshift:
    column: Z_CLUSTER
    bins: [0.2, 1.6]  # Single bin

# Example: All galaxies only (no type split)
galaxy_bins: {}  # Empty = use all galaxies
```

Then regenerate jobs:
```bash
python scripts/generate_jobs.py --config configs/rr2_clustering.yaml --array
sbatch jobs/submit_array.sh
```

## Next Steps

After results are generated:

1. **Load and analyze results programmatically**:
   ```python
   import pickle
   import numpy as np
   from pathlib import Path

   results_dir = Path('outputs_rr2')
   for result_file in results_dir.glob('*.pkl'):
       with open(result_file, 'rb') as f:
           data = pickle.load(f)

       r = data['results']['r']
       xi = data['results']['xi']
       sigma = data['results']['sigma_xi']

       # Your analysis here
   ```

2. **Compare with your existing results** in `/sps/euclid/Users/cmurray/RR2_alt/scripts/results/`

3. **Stack results** across redshift or richness bins

4. **Fit halo models** to the measured profiles

## Contact

For issues with the pipeline, check:
- [README.md](README.md) - Full documentation
- [QUICKSTART.md](QUICKSTART.md) - Quick reference
- [PIPELINE_OVERVIEW.md](PIPELINE_OVERVIEW.md) - Technical details
