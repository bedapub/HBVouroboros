Build a Singularity image of the Pisces program
===
Jitao David Zhang, July 16, 2020

There are multiple ways to install Pisces (documented in [../HBVouroboros/pisces/READMD.md](../HBVouroboros/pisces/READMD.md). Here we use a Singularity image using the following recipe.

```bash
## On HPCs using EasyBuild, we load the Singularity 3.5 with the following command
ml purge
ml load Singularity/3.5.0

singularity build pisces_s3.5.simg docker://astewart/pisces
## test it
singularity shell pisces_s3.5.simg
dotnet /app/Pisces_5.2.9.122/Pisces.dll
```

## Licensing information

The Pisces software is used with the GPL-v3 license.
