FROM python:3.12.6-slim-bullseye
MAINTAINER "Irenaeus Chan <chani@wustl.edu>"

# Volumes
VOLUME /build

ARG TAG=v0.0.1
ARG PYTHON_VERSION=3.12.6

# Install system dependencies
RUN apt-get update -qq \
    && apt-get -y install \
    build-essential \
    git \
    less \
    procps \
    libnss-sss \
    libcurl4-openssl-dev \
    curl \
    wget \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    pkg-config \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    zlib1g-dev \
    --no-install-recommends \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install packages
RUN pip install --upgrade pip \
    && pip install "setuptools<70" wheel

# Install scientific stack
RUN pip install --no-cache-dir \
    pytz \
    tzdata \
    python-dateutil \
    numpy \
    pandas \
    pysam \
    click \
    pyranges

RUN pip install --upgrade git+https://github.com/IrenaeusChan/DNA-RECAP.git
