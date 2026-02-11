FROM debian:bullseye-slim
MAINTAINER "Irenaeus Chan <chani@wustl.edu>"

# Volumes
VOLUME /build

ARG TAG=v0.0.0
ARG PYTHON_VERSION=3.12.6

# Install build dependencies
RUN apt-get update -qq \
    && apt-get -y install apt-transport-https ca-certificates \
    && apt-get update -qq \
    && apt-get -y install \
    build-essential \
    git \
    less \
    procps \
    libnss-sss \
    libcurl4-openssl-dev \
    curl \
    wget \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    libffi-dev \
    --no-install-recommends

RUN apt-get clean all

# Download and install Python
RUN wget https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz \
    && tar -xzf Python-$PYTHON_VERSION.tgz \
    && cd Python-$PYTHON_VERSION \
    && ./configure --prefix=/opt/python-$PYTHON_VERSION \
    && make \
    && make install

# Update PATH
ENV PATH="/opt/python-$PYTHON_VERSION:$PATH"

# Install pip and the package
RUN /opt/python-$PYTHON_VERSION/bin/pip3 install --upgrade pip
RUN /opt/python-$PYTHON_VERSION/bin/pip3 install --upgrade git+https://github.com/IrenaeusChan/DNA-RECAP.git
