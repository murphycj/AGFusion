"""
Download and build the agfusion database for a range of releases.
"""

import os
import subprocess

import boto3
from botocore.exceptions import ClientError

S3_BUCKET = "agfusion"
PFAM_FILE = "Pfam-A.clans.tsv"

s3_client = boto3.client("s3")


def s3_key_exists(bucket, key):
    """Return True if the given S3 key exists in the bucket."""
    try:
        s3_client.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError:
        return False


def s3_upload(local_path, bucket, key):
    """Upload a local file to S3."""
    print(f"Uploading {local_path} to s3://{bucket}/{key}")
    s3_client.upload_file(local_path, bucket, key)


def s3_download(bucket, key, local_path):
    """Download a file from S3 to a local path."""
    print(f"Downloading s3://{bucket}/{key} to {local_path}")
    s3_client.download_file(bucket, key, local_path)


def run_command(command):
    """Run a command and print it."""
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True)


# Ensure the Pfam clans file is available locally
if not os.path.exists(PFAM_FILE):
    s3_download(S3_BUCKET, PFAM_FILE, PFAM_FILE)

# species = "homo_sapiens"
species = "mus_musculus"

# Loop over releases from 96 to 110
for i in range(92, 112):
    db_key = f"agfusion.{species}.{i}.db.gz"

    # Check if the file exists on S3
    if s3_key_exists(S3_BUCKET, db_key):
        print(f"File for release {i} already exists in S3. Skipping...")
        continue

    print(f"Building release {i}")

    # Install the release using pyensembl
    run_command(f"pyensembl install --release {i} --species {species}")

    # Build the agfusion database
    run_command(f"agfusion build -d . -s {species} -r {i} --pfam {PFAM_FILE}")

    # Compress the database
    run_command(f"gzip agfusion.{species}.{i}.db")

    # Upload the compressed file to S3
    s3_upload(db_key, S3_BUCKET, db_key)

    # Delete all files for the release from pyensembl
    run_command(f"pyensembl delete-all-files --release {i} --species {species}")
