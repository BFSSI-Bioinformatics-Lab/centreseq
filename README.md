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

### Issues
The conda distribution might have an issue installing the correct version of `minced`,
 a Prokka dependency, so you might need to run the following command within your conda environment after installing centreseq:

`conda install minced`

### External Dependencies

These programs will be automatically installed with the conda package.

- [Prokka](https://github.com/tseemann/prokka)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [SeqKit](https://github.com/shenwei356/seqkit)
- [MUSCLE](https://www.drive5.com/muscle/)
- [SNP-sites](https://github.com/sanger-pathogens/snp-sites)
- [cyvcf2](https://github.com/brentp/cyvcf2)