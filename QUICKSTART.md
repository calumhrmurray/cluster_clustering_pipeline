# Quick Start Guide

This guide will get you running clustering analysis in 5 minutes.

## Step 1: Create Configuration (1 min)

```bash
cd cluster_clustering_pipeline
python scripts/create_config.py --output my_config.yaml
```

Edit `my_config.yaml` and update these paths:

```yaml
catalogues:
  clusters: /path/to/your/cluster_catalogue.fits
  galaxies: /path/to/your/galaxy_catalogue.fits
  randoms: /path/to/your/random_catalogue.fits
```

**Tip**: Use the example configs in `configs/` as templates for common use cases.

## Step 2: Test the Pipeline (1 min)

Run a single test job to verify everything works:

```bash
python scripts/run_analysis.py --config my_config.yaml --test
```

This will:
- Load your catalogues
- Run the first cluster-galaxy bin combination
- Save results to `outputs/`

Check the output file is created successfully.

## Step 3: Generate SLURM Jobs (1 min)

Create batch job scripts for all bin combinations:

```bash
python scripts/generate_jobs.py --config my_config.yaml --array
```

This creates:
- `jobs/specs/` - Job specifications (one per bin combination)
- `jobs/submit_array.sh` - SLURM submission script
- `jobs/logs/` - Directory for job logs

**Optional**: Review the generated job specs to see what will run:
```bash
ls jobs/specs/
cat jobs/specs/job_0_spec.yaml
```

## Step 4: Submit Jobs (30 sec)

Submit all jobs as a SLURM array:

```bash
sbatch jobs/submit_array.sh
```

Monitor job progress:
```bash
squeue -u $USER
```

Check logs:
```bash
tail -f jobs/logs/array_*_0.out  # First job
ls jobs/logs/                     # All logs
```

## Step 5: Plot Results (1 min)

Once jobs complete, plot the results:

```bash
python scripts/plot_results.py outputs/*.pkl --output clustering_results.png
```

For comparison plots:
```bash
python scripts/plot_results.py \
  outputs/clustering_c_*_0.20_0.50*.pkl \
  outputs/clustering_c_*_0.50_0.80*.pkl \
  --compare --output redshift_comparison.png
```

## Common Configurations

### Galaxy Type Comparison

Want to compare early-type vs late-type galaxy clustering?

```yaml
galaxy_bins:
  galaxy_type:
    column: sersic_sersic_vis_index
    bins: [0.0, 2.0, 10.0]  # Split at Sersic index = 2
```

### Multiple Redshift Bins

```yaml
cluster_bins:
  redshift:
    bins: [0.2, 0.5, 0.8, 1.2]  # 3 bins

galaxy_bins:
  redshift:
    bins: [0.0, 0.5, 1.0, 1.5, 2.0]  # 4 bins

# Total jobs = 3 × 4 = 12
```

### Weak Lensing Setup

For tangential shear measurements:

```yaml
analysis_parameters:
  mode: 2d               # Angular coordinates
  metric: Euclidean
  sep_units: deg
  min_sep: 0.01         # degrees
  max_sep: 0.5

# Lenses at low-z, sources at high-z
cluster_bins:
  redshift:
    bins: [0.2, 0.5]

galaxy_bins:
  redshift:
    bins: [0.8, 2.0]     # Must be > cluster z
```

## Troubleshooting

**Jobs fail immediately?**
- Check SLURM logs in `jobs/logs/`
- Verify catalogue paths exist
- Test locally first with `--test` flag

**No objects after filtering?**
- Run with `--dry-run` to see filter specs
- Check your bin ranges match your data
- Use `--test` mode to debug

**Out of memory?**
- Increase `mem` in SLURM config
- Reduce catalogue size for testing
- Use more CPUs if available

**Results look wrong?**
- Check coordinate columns are correct in config
- Verify cosmology parameters
- Plot with `--linear` to check non-log scale

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Explore example configs in `configs/`
- Add systematic weights for more accurate results
- Combine results from multiple bins for final analysis

## Getting Help

Common issues:
1. **"Module not found"** → Install requirements: `pip install -r requirements.txt`
2. **"File not found"** → Use absolute paths in config
3. **"TreeCorr warning"** → Usually safe to ignore patch warnings
4. **Empty output** → Check filters aren't too restrictive

Still stuck? Check the Troubleshooting section in README.md
