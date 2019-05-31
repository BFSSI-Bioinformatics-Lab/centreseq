### External Dependencies

The following must be available in your *$PATH*
- [Prokka](https://github.com/tseemann/prokka)
- [mmseqs2](https://github.com/soedinglab/MMseqs2)
- [muscle](https://www.drive5.com/muscle/)

### Installation
1. Clone the repo and navigate to directory

```
git clone https://github.com/bfssi-forest-dussault/centreseq.git
cd centreseq
```

2. Create a new conda env with the environment.yml file

```
conda env create -f environment.yml
```

3. Call the script
```
python core_pipeline.py --help
```

### Usage
```
Usage: centreseq.py [OPTIONS]

  Given an input directory containing assemblies, this will establish a core
  genome through running Prokka and mmseqs2.

Options:
  -f, --fasta-dir PATH    Path to directory containing *.fasta files to run
                          core pipeline on  [required]
  -o, --outdir PATH       Root directory to store all output files  [required]
  -n, --n-cpu INTEGER     Number of CPUs to dedicate to parallelizable steps
                          of the pipeline. Will take all available CPUs - 1 if
                          not specified.
  -m, --min-seq-id FLOAT  Sets the mmseqs cluster parameter "--min-seq-id".
                          Defaults to 0.95.
  --help                  Show this message and exit.
```