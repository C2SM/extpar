# Use the latest Ubuntu image
FROM ubuntu:latest

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update the package list and install required packages
RUN apt-get update && \
    apt-get install -y \
    git \
    vim \
    gcc \
    gfortran \
    libnetcdf-dev \
    libnetcdff-dev \
    libsz2 \
    libomp-dev \
    python3 \
    python3-pip \
    python3-venv  \
    bc \
    cdo \
    wget \
    && apt-get clean

# Verify installations
RUN gcc --version && \
    gfortran --version && \
    git --version && \
    vim --version && \
    nc-config --version && \
    nf-config --version

# Set the working directory
WORKDIR /workspace

# Copy the configuration script into the container
COPY . /workspace/

# Run the configuration script
RUN /workspace/configure.docker.gcc

# Build extpar
RUN cd /workspace && make -j 8

# Create a virtual environment and install the package
RUN python3 -m venv /workspace/venv && \
    . /workspace/venv/bin/activate && \
    cd /workspace && \
    pip install -r requirements.txt && \
    pip install setuptools && \
    python3 setup.py sdist && \
    pip install dist/extpar-*.tar.gz

# Add the virtual environment to the PATH
ENV PATH="/workspace/venv/bin:$PATH"

# Set the working directory
WORKDIR /workspace/python

# Default command
CMD ["/bin/bash"]
