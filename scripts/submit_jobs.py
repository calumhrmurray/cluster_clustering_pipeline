#!/usr/bin/env python3
"""
Generate and submit SLURM jobs in one step.
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import PipelineConfig


def main():
    parser = argparse.ArgumentParser(
        description='Generate and submit SLURM jobs for clustering analysis')
    parser.add_argument('--config', type=str, required=True,
                       help='Path to configuration file')
    parser.add_argument('--array', action='store_true', default=True,
                       help='Use job array (default: True)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Generate jobs but do not submit')

    args = parser.parse_args()

    print("="*70)
    print("GENERATE AND SUBMIT SLURM JOBS")
    print("="*70)

    # Get pipeline directory (go up from config file if in configs/)
    config_path = Path(args.config).resolve()
    if config_path.parent.name == 'configs':
        pipeline_dir = config_path.parent.parent
    else:
        pipeline_dir = config_path.parent

    # Load config to show summary
    config = PipelineConfig(args.config)
    job_specs = config.get_all_job_specs()

    print(f"\nConfiguration: {args.config}")
    print(f"Total jobs: {len(job_specs)}")
    print(f"Output directory: {config.get_output_dir()}")
    print(f"Job type: {'Array' if args.array else 'Individual'}")

    # Generate jobs by calling generate_jobs.py
    print("\n" + "-"*70)
    print("STEP 1: Generating job scripts...")
    print("-"*70)

    generate_cmd = [
        sys.executable,
        str(Path(__file__).parent / 'generate_jobs.py'),
        '--config', str(args.config)
    ]

    if args.array:
        generate_cmd.append('--array')

    result = subprocess.run(generate_cmd, capture_output=False)

    if result.returncode != 0:
        print("\n❌ Job generation failed!")
        return 1

    print("\n✓ Job scripts generated successfully")

    # Submit jobs
    if args.dry_run:
        print("\n" + "-"*70)
        print("DRY RUN - Jobs not submitted")
        print("-"*70)

        if args.array:
            submit_script = pipeline_dir / "jobs" / "submit_array.sh"
            print(f"\nTo submit, run:")
            print(f"  sbatch {submit_script}")
        else:
            submit_script = pipeline_dir / "jobs" / "scripts" / "submit_all.sh"
            print(f"\nTo submit, run:")
            print(f"  bash {submit_script}")

        return 0

    print("\n" + "-"*70)
    print("STEP 2: Submitting jobs to SLURM...")
    print("-"*70)

    if args.array:
        submit_script = pipeline_dir / "jobs" / "submit_array.sh"
        submit_cmd = ['sbatch', str(submit_script)]
    else:
        submit_script = pipeline_dir / "jobs" / "scripts" / "submit_all.sh"
        submit_cmd = ['bash', str(submit_script)]

    print(f"\nRunning: {' '.join(submit_cmd)}")

    result = subprocess.run(submit_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"\n❌ Job submission failed!")
        print(f"Error: {result.stderr}")
        return 1

    print(result.stdout)

    # Parse job ID from sbatch output
    if args.array and 'Submitted batch job' in result.stdout:
        job_id = result.stdout.strip().split()[-1]
        print(f"\n✓ Jobs submitted successfully!")
        print(f"\nJob array ID: {job_id}")
        print(f"\nMonitor with:")
        print(f"  squeue -u $USER")
        print(f"  squeue -j {job_id}")
        print(f"\nCheck logs:")
        print(f"  tail -f jobs/logs/array_{job_id}_0.out  # First job")
        print(f"  ls jobs/logs/array_{job_id}_*.out       # All jobs")
    else:
        print(f"\n✓ Jobs submitted!")
        print(f"\nMonitor with:")
        print(f"  squeue -u $USER")

    print(f"\nResults will be saved to:")
    print(f"  {config.get_output_dir()}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
