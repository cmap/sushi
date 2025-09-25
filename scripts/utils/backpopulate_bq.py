import argparse
import os
import subprocess
import sys

# The script paths should be relative to the 'scripts' directory inside your workspace
SCRIPTS_TO_RUN = ["qc_tables/qc_tables.sh", "load_bq/load_bq.sh"]


def run_backfill(build_dirs_file, workspace_path):
    """
    Reads a file containing build directories and runs specific analysis scripts for each.
    """
    print(f"🚀 Starting backfill process...")

    # --- 1. Validate paths ---
    if not os.path.isfile(build_dirs_file):
        print(f"❌ Error: The file '{build_dirs_file}' was not found.")
        sys.exit(1)

    launch_script_path = os.path.join(workspace_path, "scripts", "launch_job.sh")
    if not os.path.isfile(launch_script_path):
        print(
            f"❌ Error: The launch script 'launch_job.sh' was not found in '{os.path.join(workspace_path, 'scripts')}'"
        )
        sys.exit(1)

    # --- 2. Read the list of directories ---
    try:
        with open(build_dirs_file, "r") as f:
            build_dirs = [line.strip() for line in f if line.strip()]
    except IOError as e:
        print(f"❌ Error reading file '{build_dirs_file}': {e}")
        sys.exit(1)

    if not build_dirs:
        print("🟡 Warning: The provided file is empty. No directories to process.")
        return

    print(f"Found {len(build_dirs)} directories to process.\n")

    # --- 3. Iterate and process each directory ---
    for i, build_dir in enumerate(build_dirs, 1):
        print(f"--- ({i}/{len(build_dirs)}) Processing Directory: {build_dir} ---")

        if not os.path.isdir(build_dir):
            print(f"  🟡 Skipping: Directory does not exist.\n")
            continue

        if not os.path.isfile(os.path.join(build_dir, "config.json")):
            print(f"  🟡 Skipping: 'config.json' not found in this directory.\n")
            continue

        # --- 4. Run the designated scripts ---
        for script_name in SCRIPTS_TO_RUN:
            print(f"  ▶️  Running script: {script_name}")

            command = [launch_script_path, script_name]

            try:
                # --- THIS IS THE FIX ---
                # Create a copy of the current environment
                proc_env = os.environ.copy()

                # Set the BUILD_DIR variable for the data path
                proc_env["BUILD_DIR"] = os.path.abspath(build_dir)

                # Set the WORKSPACE variable for the code path
                proc_env["WORKSPACE"] = os.path.abspath(workspace_path)

                # Pass this complete environment to the subprocess
                result = subprocess.run(
                    command, env=proc_env, check=True, capture_output=True, text=True
                )
                print(f"  ✅ Success!")

            except subprocess.CalledProcessError as e:
                print(f"  ❌ Error executing script '{script_name}' in '{build_dir}'.")
                print(f"     Return Code: {e.returncode}")
                print(f"     Stderr:\n{e.stderr.strip()}")
                break

        print(f"--- Finished processing {build_dir} ---\n")

    print("✅ Backfill process complete.")


def main():
    """Main function to parse arguments and call the backfill runner."""
    parser = argparse.ArgumentParser(
        description="Run QC and BigQuery loading scripts for a list of build directories.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "build_dirs_file",
        type=str,
        help="Path to the text file with build directories.",
    )
    parser.add_argument(
        "--workspace",
        type=str,
        required=True,
        help="Path to your pipeline code checkout.",
    )

    args = parser.parse_args()
    run_backfill(args.build_dirs_file, args.workspace)


if __name__ == "__main__":
    main()
