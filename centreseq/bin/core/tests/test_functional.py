from centreseq import centreseq_cli
from pathlib import Path
import shutil
import tempfile

DATA_DIR = centreseq_cli.ROOT_DIR / 'data'


def test_core_and_tree_modules():
    """
    Minimal functional test to ensure the core module and tree modules are behaving as expected
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        assert tmpdir.exists()

        fastas = [
            DATA_DIR / 'sample1.fasta',
            DATA_DIR / 'sample2.fasta',
            DATA_DIR / 'sample3.fasta',
        ]

        for f in fastas:
            shutil.copy(f, str(tmpdir))

        outdir = tmpdir / 'outdir'
        outdir.mkdir()

        centreseq_cli.core_pipeline(fasta_dir=tmpdir, outdir=outdir,
                                    n_cpu=centreseq_cli.set_cpu_count(), n_cpu_pickbest=2,
                                    min_seq_id=0.95, coverage_length=0.95, medoid_repseqs=False, pairwise=True)

        # Make sure output directory looks like we expect
        assert (outdir / 'reports').exists()
        assert (outdir / 'core_genome').exists()
        assert (outdir / 'mmseqs2').exists()
        assert (outdir / 'prokka').exists()

        # Pairwise output
        assert (outdir / 'static').exists()
        assert (outdir / 'network_graph.html').exists()
        assert (outdir / 'network_graph_coding.tsv').exists()

        # Check contents of the core gene count file, make sure expected #s are correct
        core_gene_report = outdir / 'reports' / 'core_gene_count_report.txt'
        assert core_gene_report.exists()
        with open(core_gene_report, 'r') as f:
            lines = list(f.readlines())
            assert lines[0].split()[-1] == '1'  # Expected =100% count
            assert lines[1].split()[-1] == '2'  # Expected >=90% count
            assert lines[2].split()[-1] == '5'  # Expected >=50% count

        # Testing the tree module with the output from the core module
        treedir = outdir / 'tree'
        treedir.mkdir()
        centreseq_cli.tree_pipeline(
            summary_report=outdir / 'reports' / 'summary_report.tsv',
            prokka_dir=outdir / 'prokka',
            outdir=treedir,
            n_cpu=centreseq_cli.set_cpu_count(),
            percentile=100.0
        )

        shutil.rmtree("/home/brock/Documents/testdir")
        shutil.copytree(str(outdir), "/home/brock/Documents/testdir")

        assert (treedir / 'aligned_loci').exists()
        assert (treedir / 'concatenated_sequences').exists()

        treefile = treedir / 'concatenated_sequences' / 'concatenated_sequences.fasta'
        assert treefile.exists()

        # Open up the report file and check if there is only 1 unique sequence represented
        with open(str(treefile), 'r') as f:
            seqs = []
            for line in f:
                if line[0] == ">":
                    continue
                seqs.append(line.strip())
            assert len(set(seqs)) == 1
