[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/centreseq/README.html)
[![Documentation Status](https://readthedocs.org/projects/centreseq/badge/?version=latest)](https://centreseq.readthedocs.io/en/latest/?badge=latest)

The most recent documentation can be found at the following URL: https://centreseq.readthedocs.io/en/latest/

### Quick installation in a Conda environment
Due to some strange Conda/Bioconda issues with both centreseq and Prokka, 
the installation process is convoluted to get everything working. The following incantation should work: 
```text
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n centreseq python=3.7
conda activate centreseq
conda install centreseq -c bioconda
conda uninstall prokka
conda install -c conda-forge -c bioconda -c defaults prokka
conda install centreseq -c bioconda
```


### Developer notes/Known issues
- Version of mmseqs2 has been pinned at 9-d36de in meta.yaml due to 10-6d92c containing a bug which causes the functional test to 
hang indefinitely on some systems.
- Prokka's installation appears to be broken in some cases when just running `conda install prokka -c bioconda`, 
the solution is to downgrade Perl with `conda install perl=5.22.0`