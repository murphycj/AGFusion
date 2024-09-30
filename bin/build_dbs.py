"""
Download and build the agfusion database for a range of releases.
"""

import subprocess


def run_command(command):
    """Run a command and print it."""
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)


species = "homo_sapiens"

# Loop over releases from 96 to 110
for i in range(96, 112):
    # Check if the file exists on S3
    s3_check_command = f"aws s3 ls s3://agfusion/agfusion.{species}.{i}.db.gz"
    result = subprocess.run(
        s3_check_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )

    if result.returncode == 0:
        print(f"File for release {i} already exists in S3. Skipping...")
        continue

    # continue

    print(f"Building release {i}")

    # Install the release using pyensembl
    run_command(f"pyensembl install --release {i} --species {species}")

    # Build the agfusion database
    run_command(f"agfusion build -d . -s {species} -r {i} --pfam Pfam-A.clans.tsv")

    # Compress the database
    run_command(f"gzip agfusion.{species}.{i}.db")

    # Upload the compressed file to S3
    run_command(f"aws s3 cp agfusion.{species}.{i}.db.gz s3://agfusion")

    # Delete all files for the release from pyensembl
    run_command(f"pyensembl delete-all-files --release {i} --species {species}")
