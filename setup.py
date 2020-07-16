from setuptools import setup, find_packages

setup(
    name='HBVouroboros',
    version='0.6.4',
    description="sequence-based HBV genotyping and expression profiling",
    packages=find_packages(),

    license="GPL-3",
    keywords="HBV bioinformatics sequencing genotyping expression",
    # the following three lines are not required
    author = 'Jitao David Zhang',
    url = 'https://github.roche.com/BEDA/HBVouroboros',
    # python modules
    py_modules = ['HBVouroboros'],
    # names of binaries
    scripts = [
        'bin/HBVouroboros_map_samples.py',
        'bin/HBVouroboros_map_biokit.py',
        'bin/HBVouroboros_build_refgenomes.py',
        'bin/HBVouroboros_trimmomatic.py'
        ],
    package_data = {
        "HBVouroboros": ["build_refgenomes/Snakefile", 
            "config/cluster.json",
            "align_reads/Snakefile",
            "align_reads/map-fastq-snakefile",
            "align_reads/config/config.yaml",
            "trim_reads/Snakefile",
            "trim_reads/config/primers.fasta",
            "trim_reads/config/config.yaml",
            "biokit/Snakefile"]
    }
    # install_requires = ['requests>=2.19.1']
)

