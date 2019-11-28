[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/centreseq/README.html)
[![Documentation Status](https://readthedocs.org/projects/centreseq/badge/?version=latest)](https://centreseq.readthedocs.io/en/latest/?badge=latest)

The most recent documentation can be found at the following URL: https://centreseq.readthedocs.io/en/latest/

### Quick installation in a Conda environment
```text
conda create -n centreseq python=3.7
conda activate centreseq
conda install centreseq -c bioconda
```

### Developer notes
- Version of mmseqs2 has been pinned at 9-d36de in meta.yaml due to 10-6d92c containing a bug which causes the functional test to 
hang indefinitely on some systems.