# Define base image
FROM continuumio/miniconda3
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

# Set working directory for the project
WORKDIR /app

# Create Conda environment from YAML file
RUN echo "Cloning HBVouroboros and setting up environment..."
RUN git clone -b final_corrections https://github.com/bedapub/HBVouroboros.git
RUN conda env create -f ./HBVouroboros/envs/environment.yml
SHELL ["conda", "run", "-n", "HBVouroboros", "/bin/bash", "-c"]
