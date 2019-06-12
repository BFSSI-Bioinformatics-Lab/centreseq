[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/centreseq/README.html)


### Installation
1. Create a new conda environment (*optional, but strongly recommended*)

```
conda env create -n centreseq python=3.6
conda activate centreseq
```

2. Install via bioconda
```
conda install -c bioconda centreseq
```

OR

Install with docker
```
docker pull quay.io/biocontainers/centreseq:<tag>
```

### Usage
```
Usage: centreseq [OPTIONS] COMMAND [ARGS]...

  centreseq builds an annotated core genome using assemblies as input.

Options:
  --version   Print the version and exit.
  --citation  Print the citation for this software and exit.
  --help      Show this message and exit.

Commands:
  core  Given an input directory containing assemblies, establishes a core
        genome
  tree  Produces output for phylogenetic tree software
```

### External Dependencies

These programs will be automatically installed with the conda package.

- [Prokka](https://github.com/tseemann/prokka)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [SeqKit](https://github.com/shenwei356/seqkit)
- [MUSCLE](https://www.drive5.com/muscle/)
- [SNP-sites](https://github.com/sanger-pathogens/snp-sites)
- [cyvcf2](https://github.com/brentp/cyvcf2)