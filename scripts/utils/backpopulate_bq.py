import argparse
import os
import subprocess
import json
import sys

# Define the specific scripts you want to run for the backfill.
# These paths are relative to the 'scripts' directory in your workspace.
SCRIPTS_TO_RUN = ["../qc_tables/qc_tables.sh", "../load_bq/load_bq.sh"]


def run_backfill(build_dirs_file, workspace_path):
    """
    Reads a file containing build directories and runs specific analysis scripts for each.

    Args:
        build_dirs_file (str): Path to a text file with one build directory per line.
        workspace_path (str): Path to the root of your pipeline code checkout
                              (the directory containing the 'scripts' folder).
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
        print(
            "Please ensure the --workspace argument points to the correct pipeline repository."
        )
        sys.exit(1)

    # --- 2. Read the list of directories ---
    try:
        with open(build_dirs_file, "r") as f:
            # Read all lines, stripping whitespace and filtering out empty lines
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

        # Validate that the build directory and its config.json exist
        if not os.path.isdir(build_dir):
            print(f"  🟡 Skipping: Directory does not exist.\n")
            continue

        config_path = os.path.join(build_dir, "config.json")
        if not os.path.isfile(config_path):
            print(f"  🟡 Skipping: 'config.json' not found in this directory.\n")
            continue

        # The scripts rely on the config, but we don't need to read it in Python.
        # The launch_job.sh script will handle it, provided we run from the correct directory.

        # --- 4. Run the designated scripts ---
        for script_name in SCRIPTS_TO_RUN:
            print(f"  ▶️  Running script: {script_name}")

            command = [launch_script_path, script_name]

            try:
                # The key is to execute the command with the 'current working directory' (cwd)
                # set to the build_dir. This allows the scripts to find 'config.json'.
                result = subprocess.run(
                    command,
                    cwd=build_dir,
                    check=True,  # Raises an exception if the command returns a non-zero exit code
                    capture_output=True,  # Captures stdout and stderr
                    text=True,  # Decodes stdout/stderr as text
                )
                print(f"  ✅ Success!")
                # Optional: print stdout if you want to see the script's output
                # if result.stdout:
                #     print(f"     Output: {result.stdout.strip()}")

            except FileNotFoundError:
                print(
                    f"  ❌ Error: Command not found. Ensure '{launch_script_path}' is a valid executable."
                )
                break  # Stop processing this directory
            except subprocess.CalledProcessError as e:
                print(f"  ❌ Error executing script '{script_name}' in '{build_dir}'.")
                print(f"     Return Code: {e.returncode}")
                print(f"     Stderr:\n{e.stderr.strip()}")
                break  # Stop processing this directory and move to the next one
            except Exception as e:
                print(f"  ❌ An unexpected error occurred: {e}")
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
        help="Path to the text file containing the list of build directories, one per line.",
    )

    parser.add_argument(
        "--workspace",
        type=str,
        required=True,
        help="Path to your local checkout of the pipeline code.\n"
        "This directory should contain the 'scripts/launch_job.sh' file.",
    )

    args = parser.parse_args()
    run_backfill(args.build_dirs_file, args.workspace)


if __name__ == "__main__":
    main()
