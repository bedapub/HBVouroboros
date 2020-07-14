Detecting SNVs with Pisces and HBVouroboros
=====

## Background

To use [Illumina/Pisces](https://github.com/Illumina/Pisces/) for SNV calling, there are two options

1. [Use the Bioconda package](#working-with-the-pisces-package-in-bioconda)
   * Advantage: no need to work with Docker images and to mount data
   * Disadvantage: need to install dotnet runtime 2.x
1. [Build and use a Docker image](#working-with-a-docker-image-of-pisces)
   * Advantage: no need to install dotnet runtime
   * Disadvantage: Both `sudo` permission and Docker installation are required.

Below we document both possibilities. During the development of *HBVouroboros*,
I used the Bioconda option. In case you identify issues when using the Docker
option, please kindly let me know.

## Working with the Pisces package in Bioconda

The package [`Pisces`](https://anaconda.org/bioconda/pisces) is already
installed as part of the `HBVouroboros` environment of conda. The latest version
as of July 2020 is 5.2.10.49.

To use `Pisces`, activate the conda environment (`conda activate HBVouroboros`)
and run `pisces` on the command line. If the `dotnet` runtime is not installed,
you will receive a complaint.

You have to download and install the .NET (dotnet) runtime before running
`Pisces` from [dotnet.microsoft.com](https://dotnet.microsoft.com/download).
Important: though the latest dotnet runtime is of version 3.1, `Pisces` requires
version 2.x (for instance 2.1).

Depending on the operation system you are using, the command to install dotnet
varies. In my setting (64-bit Mint Linux), I used following commands (which I
discovered from a reply on
[StackOverflow](https://stackoverflow.com/questions/52737293/install-dotnet-core-on-linux-mint-19):

```bash
wget https://packages.microsoft.com/config/debian/10/packages-microsoft-prod.deb \
    -O packages-microsoft-prod.deb
sudo dpkg -i packages-microsoft-prod.deb
sudo apt-get update
sudo apt-get install -y apt-transport-https
sudo apt-get update
sudo apt-get install -y dotnet-runtime-2.1 ## Pisces does not work yet with 3.1
```

## Working with a Docker image of Pisces

Two options are here: [building your own Docker
image](#building-the-image-from-source-code), or [using the image available on
Docker Hub](#working-with-pre-built-Docker-image).

### Building the image from source code

The Dockerfile and run_analysis.sh is copied from the [GitHub repository of
Pisces](https://github.com/Illumina/Pisces/tree/master/docker/ExamplePiscesPaperAnalysis),
sub-folder `docker/ExamplePisces5.2.5Release`, version 5.2.5, commit
[5dce4fe](https://github.com/Illumina/Pisces/commit/5dce4fe7d1dc4603ca35affe258cbce14cf4ae1c).

## Instructions

To build the Docker image, run the following command in the `config` directory.

```bash
sudo docker build -t pisces . > docker_build.log
```

To run the command in the simple mode, run

```bash
sudo docker run pisces
```

To run Pisces with mounted data, we need to first mount the data, and inspect
with bash

```bash
docker run -it --mount type=bind,source="$(pwd)"/data,target=/data pisces
/bin/bash
```

### Working with pre-built Docker image

An Docker image of pisces is available on [Docker
Hub](https://hub.docker.com/r/astewart/pisces). It can be pulled by Docker by

```bash
docker pull astewart/pisces
```

If Singularity is preferred, use the following commands to convert the Docker
image into a Singularity one, and run the Pisces binary program from the
Singularity image within.

```bash
singularity build pisces.simg docker://astewart/pisces
singularity shell pisces.simg
dotnet /app/Pisces_5.2.9.122/Pisces.dll
```

## Licensing information

The Pisces software is used with the GPL-v3 license.
