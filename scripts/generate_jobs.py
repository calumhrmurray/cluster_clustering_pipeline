#!/usr/bin/env python3
"""
Generate SLURM batch job scripts for parallel clustering analysis.
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path to import src modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import PipelineConfig


def create_job_script(job_spec, config, script_dir, pipeline_dir):
    """
    Create a SLURM batch script for a single job.

    Parameters
    ----------
    job_spec : dict
        Job specification with cluster and galaxy filters
    config : PipelineConfig
        Pipeline configuration
    script_dir : Path
        Directory to save job scripts
    pipeline_dir : Path
        Path to pipeline directory

    Returns
    -------
    Path
        Path to created job script
    """
    job_id = job_spec['job_id']
    slurm_config = config.get_slurm_config()

    # Create job name
    from src.filters import CatalogueFilter
    cluster_label = CatalogueFilter.get_bin_label(job_spec['cluster_filters'])
    galaxy_label = CatalogueFilter.get_bin_label(job_spec['galaxy_filters'])
    job_name = f"cluster_clust_{job_id}"

    # SLURM script header
    script_content = "#!/bin/bash\n"
    script_content += f"#SBATCH --job-name={job_name}\n"
    script_content += f"#SBATCH --output={pipeline_dir}/jobs/logs/job_{job_id}_%j.out\n"
    script_content += f"#SBATCH --error={pipeline_dir}/jobs/logs/job_{job_id}_%j.err\n"
    script_content += f"#SBATCH --time={slurm_config.get('time', '02:00:00')}\n"
    script_content += f"#SBATCH --mem={slurm_config.get('mem', '16G')}\n"
    script_content += f"#SBATCH --cpus-per-task={slurm_config.get('cpus', 4)}\n"

    if 'partition' in slurm_config and slurm_config['partition']:
        script_content += f"#SBATCH --partition={slurm_config['partition']}\n"

    if 'account' in slurm_config and slurm_config['account']:
        script_content += f"#SBATCH --account={slurm_config['account']}\n"

    # Additional SLURM options
    for key, value in slurm_config.items():
        if key not in ['time', 'mem', 'cpus', 'partition', 'account'] and value:
            script_content += f"#SBATCH --{key}={value}\n"

    script_content += "\n"

    # Environment setup
    script_content += "# Environment setup\n"
    script_content += "set -e  # Exit on error\n"
    script_content += "set -u  # Exit on undefined variable\n"
    script_content += "\n"

    # Load modules if specified
    if 'modules' in slurm_config and slurm_config['modules']:
        script_content += "# Load required modules\n"
        for module in slurm_config['modules']:
            script_content += f"module load {module}\n"
        script_content += "\n"

    # Activate conda/virtual environment if specified
    if 'conda_env' in slurm_config and slurm_config['conda_env']:
        script_content += "# Activate conda environment\n"
        script_content += f"source $(conda info --base)/etc/profile.d/conda.sh\n"
        script_content += f"conda activate {slurm_config['conda_env']}\n"
        script_content += "\n"
    elif 'venv' in slurm_config and slurm_config['venv']:
        script_content += "# Activate virtual environment\n"
        script_content += f"source {slurm_config['venv']}/bin/activate\n"
        script_content += "\n"

    # Change to pipeline directory
    script_content += f"cd {pipeline_dir}\n"
    script_content += "\n"

    # Run the analysis
    script_content += "# Run clustering analysis\n"
    script_content += "echo 'Starting clustering analysis...'\n"
    script_content += f"echo 'Job ID: {job_id}'\n"
    script_content += f"echo 'Cluster bin: {cluster_label}'\n"
    script_content += f"echo 'Galaxy bin: {galaxy_label}'\n"
    script_content += "echo ''\n"
    script_content += "\n"

    # Create job specification file
    job_spec_file = pipeline_dir / f"jobs/specs/job_{job_id}_spec.yaml"
    script_content += f"# Job specification is in {job_spec_file}\n"

    # Run the Python script
    config_file = pipeline_dir / "config.yaml"
    script_content += f"python scripts/run_single_job.py \\\n"
    script_content += f"    --config {config_file} \\\n"
    script_content += f"    --job-spec {job_spec_file}\n"
    script_content += "\n"

    script_content += "echo 'Analysis complete!'\n"

    # Write script to file
    script_path = script_dir / f"job_{job_id}.sh"
    with open(script_path, 'w') as f:
        f.write(script_content)

    # Make executable
    script_path.chmod(0o755)

    return script_path


def create_job_array_script(config, pipeline_dir, n_jobs, config_path):
    """
    Create a single SLURM job array script for all jobs.

    Parameters
    ----------
    config : PipelineConfig
        Pipeline configuration
    pipeline_dir : Path
        Pipeline directory
    n_jobs : int
        Total number of jobs
    config_path : Path
        Path to config file

    Returns
    -------
    Path
        Path to job array script
    """
    slurm_config = config.get_slurm_config()

    script_content = "#!/bin/bash\n"
    script_content += f"#SBATCH --job-name=cluster_clust_array\n"
    script_content += f"#SBATCH --output={pipeline_dir}/jobs/logs/array_%A_%a.out\n"
    script_content += f"#SBATCH --error={pipeline_dir}/jobs/logs/array_%A_%a.err\n"
    script_content += f"#SBATCH --array=0-{n_jobs-1}\n"
    script_content += f"#SBATCH --time={slurm_config.get('time', '02:00:00')}\n"
    script_content += f"#SBATCH --mem={slurm_config.get('mem', '16G')}\n"
    script_content += f"#SBATCH --cpus-per-task={slurm_config.get('cpus', 4)}\n"

    if 'partition' in slurm_config and slurm_config['partition']:
        script_content += f"#SBATCH --partition={slurm_config['partition']}\n"

    if 'account' in slurm_config and slurm_config['account']:
        script_content += f"#SBATCH --account={slurm_config['account']}\n"

    script_content += "\n"

    # Environment setup
    script_content += "set -e\n"
    script_content += "set -u\n"
    script_content += "\n"

    # Load modules
    if 'modules' in slurm_config and slurm_config['modules']:
        for module in slurm_config['modules']:
            script_content += f"module load {module}\n"
        script_content += "\n"

    # Activate environment
    if 'conda_env' in slurm_config and slurm_config['conda_env']:
        script_content += f"source $(conda info --base)/etc/profile.d/conda.sh\n"
        script_content += f"conda activate {slurm_config['conda_env']}\n"
        script_content += "\n"

    script_content += f"cd {pipeline_dir}\n"
    script_content += "\n"

    # Use SLURM_ARRAY_TASK_ID to select job
    script_content += "JOB_ID=$SLURM_ARRAY_TASK_ID\n"
    script_content += "echo \"Running job $JOB_ID\"\n"
    script_content += "\n"

    # Use the actual config file path
    script_content += f"python scripts/run_single_job.py \\\n"
    script_content += f"    --config {config_path} \\\n"
    script_content += f"    --job-spec jobs/specs/job_${{JOB_ID}}_spec.yaml\n"

    # Write script
    script_path = pipeline_dir / "jobs/submit_array.sh"
    with open(script_path, 'w') as f:
        f.write(script_content)

    script_path.chmod(0o755)

    return script_path


def main():
    parser = argparse.ArgumentParser(
        description='Generate SLURM job scripts for clustering analysis')
    parser.add_argument('--config', type=str, required=True,
                       help='Path to configuration file')
    parser.add_argument('--array', action='store_true',
                       help='Generate a single job array script instead of individual scripts')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Directory for job scripts (default: pipeline_dir/jobs)')

    args = parser.parse_args()

    # Load configuration
    config = PipelineConfig(args.config)

    # Get pipeline directory (go up from config file if in configs/)
    config_path = Path(args.config).resolve()
    if config_path.parent.name == 'configs':
        pipeline_dir = config_path.parent.parent
    else:
        pipeline_dir = config_path.parent

    # Create output directories
    if args.output_dir:
        job_dir = Path(args.output_dir)
    else:
        job_dir = pipeline_dir / "jobs"

    script_dir = job_dir / "scripts"
    spec_dir = job_dir / "specs"
    log_dir = job_dir / "logs"

    for directory in [script_dir, spec_dir, log_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    # Generate job specifications
    print("Generating job specifications...")
    job_specs = config.get_all_job_specs()

    # Save individual job specs to YAML files
    import yaml

    def convert_tuples_to_lists(obj):
        """Recursively convert tuples to lists for YAML compatibility."""
        if isinstance(obj, dict):
            return {k: convert_tuples_to_lists(v) for k, v in obj.items()}
        elif isinstance(obj, tuple):
            return list(obj)
        elif isinstance(obj, list):
            return [convert_tuples_to_lists(item) for item in obj]
        else:
            return obj

    for job_spec in job_specs:
        job_id = job_spec['job_id']
        spec_file = spec_dir / f"job_{job_id}_spec.yaml"

        # Convert tuples to lists for safe YAML serialization
        safe_spec = convert_tuples_to_lists(job_spec)

        with open(spec_file, 'w') as f:
            yaml.dump(safe_spec, f, default_flow_style=False)

    print(f"Created {len(job_specs)} job specifications")

    if args.array:
        # Create single job array script
        print("\nGenerating job array script...")
        array_script = create_job_array_script(config, pipeline_dir, len(job_specs), config_path)
        print(f"Job array script created: {array_script}")
        print(f"\nTo submit: sbatch {array_script}")

    else:
        # Create individual job scripts
        print("\nGenerating individual job scripts...")
        for job_spec in job_specs:
            script_path = create_job_script(job_spec, config, script_dir, pipeline_dir)

        print(f"Created {len(job_specs)} job scripts in {script_dir}")

        # Create submission script
        submit_all = script_dir / "submit_all.sh"
        with open(submit_all, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# Submit all jobs\n\n")
            for job_spec in job_specs:
                job_id = job_spec['job_id']
                f.write(f"sbatch {script_dir}/job_{job_id}.sh\n")

        submit_all.chmod(0o755)

        print(f"\nSubmission script created: {submit_all}")
        print(f"To submit all jobs: bash {submit_all}")

    print("\nDone!")


if __name__ == '__main__':
    main()
