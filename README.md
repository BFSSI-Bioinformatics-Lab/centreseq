### External Dependencies

The following must be available in your *$PATH*.
These programs will be automatically installed with the conda package.

- [Prokka](https://github.com/tseemann/prokka)
- [mmseqs2](https://github.com/soedinglab/MMseqs2)
- [seqkit](https://github.com/shenwei356/seqkit)
- [muscle](https://www.drive5.com/muscle/)

### Installation
1. Create a new conda environment (*optional, but strongly recommended*)

```
conda env create -n centreseq python=3.6
```

2. Install via bioconda
```
conda install -c bioconda centreseq  # Coming soon
```

### Usage
```
Usage: centreseq [OPTIONS] COMMAND [ARGS]...

  centreseq builds an annotated core genome using assemblies as input.

Options:
  --help  Show this message and exit.

Commands:
  core  Given an input directory containing assemblies, establishes a core
        genome
  tree  Produces output for phylogenetic tree software
```