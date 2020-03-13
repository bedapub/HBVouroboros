from setuptools import setup, find_packages

setup(
    name = 'HBVouroboros',
    version = '0.1',
    description = "Map reads to HBV cccDNA genomes",
    packages=find_packages(),

    # the following three lines are not required
    author = 'Jitao David Zhang',
    url = 'https://github.roche.com/BEDA/HBVouroboros',
    # python modules
    py_modules = ['HBVouroboros'],
    # names of binaries
    # scripts = ['run_bppc.py'],
    # list all the requirements here
    # install_requires = ['requests>=2.19.1']
)

