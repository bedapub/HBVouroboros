from setuptools import setup, find_packages

setup(
    name = 'HBVouroboros',
    version = '0.6',
    description = "Map reads to HBV cccDNA genomes",
    packages=find_packages(),

    # the following three lines are not required
    author = 'Jitao David Zhang',
    url = 'https://github.roche.com/BEDA/HBVouroboros',
    # python modules
    py_modules = ['HBVouroboros'],
    # names of binaries
    scripts = [
        'bin/HBVouroboros_map_biokit.py',
        'bin/HBVouroboros_build_refgenomes.py'
        ],
    package_data = {
        "HBVouroboros": ["build_refgenomes/Snakefile", 
            "align_reads/Snakefile",
            "align_reads/config/cluster.json",
            "align_reads/config/config.yaml",
            "biokit/Snakefile"
            ]}
    # list all the requirements here
    # install_requires = ['requests>=2.19.1']
)

