FROM python:3.10-slim

LABEL maintainer="Charlie Murphy <charliemurphyj@gmail.com>" \
    description="AGFusion - annotate gene fusions"

RUN apt-get update

# Add the COMPASS source files to the container
WORKDIR /app
COPY . .

# Install COMPASS
RUN \
    pip install --quiet --upgrade pip && \
    pip install . && \
    pyensembl install --species homo_sapiens --release 75

WORKDIR /data/

# Check everything is working smoothly
RUN echo "Docker build log: Testing agfusion" 1>&2 && \
    agfusion --help \
