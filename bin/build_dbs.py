"""
Download and build the agfusion database for a range of releases.
"""

import logging
import os
import subprocess

import boto3
from botocore.exceptions import ClientError

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

S3_BUCKET = "agfusion"
PFAM_FILE = "Pfam-A.clans.tsv"

s3_client = boto3.client("s3")


def run_command(command):
    """Run a command and log it."""
    logger.info("Running: %s", command)
    subprocess.run(command, shell=True, check=True)


# Ensure the Pfam clans file is available locally
if not os.path.exists(PFAM_FILE):
    logger.info("Downloading s3://%s/%s", S3_BUCKET, PFAM_FILE)
    s3_client.download_file(S3_BUCKET, PFAM_FILE, PFAM_FILE)

# species = "homo_sapiens"
species = "mus_musculus"

# Loop over releases from 96 to 110
for i in range(92, 112):
    db_key = f"agfusion.{species}.{i}.db.gz"

    # Check if the file exists on S3
    try:
        s3_client.head_object(Bucket=S3_BUCKET, Key=db_key)
        logger.info("File for release %d already exists in S3. Skipping...", i)
        continue
    except ClientError:
        pass

    logger.info("Building release %d", i)

    # Install the release using pyensembl
    run_command(f"pyensembl install --release {i} --species {species}")

    # Build the agfusion database
    run_command(f"agfusion build -d . -s {species} -r {i} --pfam {PFAM_FILE}")

    # Compress the database
    run_command(f"gzip agfusion.{species}.{i}.db")

    # Upload the compressed file to S3
    logger.info("Uploading %s to s3://%s/%s", db_key, S3_BUCKET, db_key)
    s3_client.upload_file(db_key, S3_BUCKET, db_key)

    # Delete all files for the release from pyensembl
    run_command(f"pyensembl delete-all-files --release {i} --species {species}")
