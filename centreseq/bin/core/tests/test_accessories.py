from centreseq.bin.core.accessories import *

DATA_DIR = Path(__file__).parent / 'data'


def test_get_list_union():
    x = [1, 2, 3]
    y = [2, 3, 4]
    expected = {1, 2, 3, 4}
    result = get_list_union(x, y)
    assert expected == result


def test_extract_sample_id_from_fasta():
    input_fasta = DATA_DIR / 'tiny_genome.fasta'
    expected = 'tiny_genome'
    result = extract_sample_id_from_fasta(input_fasta)
    assert expected == result


def test_get_fasta_headers():
    input_fasta = DATA_DIR / 'tiny_genome.fasta'
    expected = ['>18-10122.pilon_00001 hypothetical protein']
    result = get_fasta_headers(input_fasta)
    assert expected == result


def test_set_cpu_count():
    expected = 4
    result = set_cpu_count(4)
    assert expected == result


def test_generated_unordered_pairs():
    values = [1, 2, 3, 4]
    expected = [(1, 2), (1, 3), (1, 4),
                (2, 3), (2, 4), (3, 4)]
    result = generate_unordered_pairs(values)
    assert expected == result


def test_sort_fasta():
    input_fasta = DATA_DIR / 'sort_test.fasta'
    sorted_fasta = DATA_DIR / 'sort_test.sorted.fasta'
    sort_fasta(input_fasta, remove_original=False)

    input_dict = SeqIO.to_dict(SeqIO.parse(str(input_fasta), "fasta"))
    sorted_dict = SeqIO.to_dict(SeqIO.parse(str(sorted_fasta), "fasta"))

    # cleanup
    sorted_fasta.unlink()

    # Ensure contents are the same after sorting procedure
    for key, val in input_dict.items():
        assert str(val) == str(sorted_dict[key])
