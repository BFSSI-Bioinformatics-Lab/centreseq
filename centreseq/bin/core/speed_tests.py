from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import gzip
import re
import tempfile
import pandas as pd
import subprocess
import csv


def faster_fasta_searching(infile, seqs_oi = []):
    """
    Loads up a fasta file into a list. Much faster than using SeqIO
    :param infile: fasta infile
    :param seqs_oi:  a list of sequence ids you want to keep. If you want to keep everything pass []
    :return:
    """
    skip = True
    gene_name = ""
    gene_description = ""
    seq = ""
    seqs_all = []
    f = open_possible_gzip(infile)
    for line in f:
        if line[0] == ">":
            # Resolve last gene
            if (seqs_oi == []) | (gene_name in seqs_oi):
                seqs_all.append(SeqRecord(Seq(seq), id=gene_name, name=gene_name, description = gene_description))
            # Initialize new gene
            seq = ""
            gene_name = line.split(" ")[0].lstrip(">")
            gene_description = line.rstrip("\n")
            if seqs_oi == []:
                skip = False
            else:
                if gene_name in seqs_oi:
                    skip = False
                else:
                    skip = True
        elif skip:
            continue
        else:
            seq += line.rstrip("\n")
    f.close()
    return seqs_all


def read_seqs(infile, filter_list=None):
    """
    Reads up sequences from a path to a fasta file
    :param infile: path to fasta file
    :param filter_list: Strings that should be in the description of the sequences
    :return: a list of strings
    """
    r = []
    f = open_possible_gzip(infile)
    for seq in SeqIO.parse(f, "fasta"):
        if filter_list is not None:
            assert isinstance(filter_list, list)
            if any([x in seq.description for x in filter_list]):
                r.append(seq)
        else:
            r.append(seq)
    f.close()
    return r


def open_possible_gzip(infile, flags="rt"):
    """
    Opens a file handle for a gzipped or non-zipped file
    :param infile: Path to file
    :param flags:
    :return: file handle
    """
    if re.search("\.gz$", infile):
        f = gzip.open(infile, flags)
    else:
        f = open(infile, flags)
    return f


def write_random_data(n_rows = 10000, n_pts = 5):
    results = []
    f_out = tempfile.NamedTemporaryFile(delete = False)
    for i in range(n_rows):
        results.append(["A" for j in range(n_pts)])
    pd.DataFrame(results).to_csv(f_out.name, sep = ",", index = False)
    return f_out.name


def read_data_csv(infile):
    r = []
    with open(infile) as f:
        reader = csv.reader(f, delimiter = ",")
        for line in reader:
            r.append(line[0])
    return r


def read_data_pandas(infile):
    return pd.read_csv(infile)


def loop_data_pandas(Df):
    r = []
    for (i,dat) in Df.iterrows():
        r.append(dat[0])
    return r


def write_seqs_to_file(seq_list, outfile_seq=None):
    """
    Write sequences to file. If not file is given then this is written to a tempfile
    :param seq_list:  a list of sequence objects
    :param outfile_seq: outfile path
    :return: the name of the output file
    """
    if outfile_seq is None:
        outfile_seq = tempfile.NamedTemporaryFile(suffix=".fasta", delete=True).name
    with open(outfile_seq, "w") as f:
        SeqIO.write(seq_list, f, "fasta")
    return outfile_seq


def run_blast(seqs1, seqs2):
    file1 = write_seqs_to_file(seqs1)
    file2 = write_seqs_to_file(seqs2)

    print(file1, file2)
    # Make a DB
    cmd = ["makeblastdb", '-in', file1, '-dbtype', 'nucl']
    subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    # BLAST
    cmd = ['blastn', "-query", file2, "-db", file1, "-outfmt", "6"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    blast_output = proc.stdout.read().decode("utf-8")
    blast_output = blast_output.rstrip("\n")
    return blast_output


def run_mmseqs(seqs1, seqs2):
    file1 = write_seqs_to_file(seqs1)
    file2 = write_seqs_to_file(seqs2)
    subprocess.call(["rm","-f","test_out.csv"])
    print(file1, file2)

    cmd = ["mmseqs", "easy-search", file1, file2, "test_out.csv", "test/"]
    subprocess.call(cmd)
    results = []
    with open("test_out.csv") as f:
        for line in csv.reader(f, delimiter = "\t"):
            results.append(line)
    return results

if __name__ == "__main__":
    cmd = ["mmseqs", "easy-search"]
    subprocess.call(cmd)

    infile = "1250Genomes_HardCorePrototype/prokka/SRS3983999.pilon/SRS3983999.pilon.ffn"

    seqs_oi = ["SRS3983999.pilon_04413" ,"SRS3983999.pilon_04414" ,"SRS3983999.pilon_04415" ,"SRS3983999.pilon_04416" ,"SRS3983999.pilon_04417" ,"SRS3983999.pilon_04418" ,"SRS3983999.pilon_04419" ,"SRS3983999.pilon_04420" ,"SRS3983999.pilon_04421" ,"SRS3983999.pilon_04422" ,"SRS3983999.pilon_04423" ,"SRS3983999.pilon_04424" ,"SRS3983999.pilon_04425" ,"SRS3983999.pilon_04426" ,"SRS3983999.pilon_04427" ,"SRS3983999.pilon_04428" ,"SRS3983999.pilon_04429" ,"SRS3983999.pilon_04430" ,"SRS3983999.pilon_04431" ,"SRS3983999.pilon_04432" ,"SRS3983999.pilon_04433" ,"SRS3983999.pilon_04434" ,"SRS3983999.pilon_04435" ,"SRS3983999.pilon_04436" ,"SRS3983999.pilon_04437" ,"SRS3983999.pilon_04438" ,"SRS3983999.pilon_04439" ,"SRS3983999.pilon_04440" ,"SRS3983999.pilon_04441" ,"SRS3983999.pilon_04442" ,"SRS3983999.pilon_04443" ,"SRS3983999.pilon_04444" ,"SRS3983999.pilon_04445" ,"SRS3983999.pilon_04446" ,"SRS3983999.pilon_04447" ,"SRS3983999.pilon_04448" ,"SRS3983999.pilon_04449" ,"SRS3983999.pilon_04450" ,"SRS3983999.pilon_04451" ,"SRS3983999.pilon_04452" ,"SRS3983999.pilon_04453" ,"SRS3983999.pilon_04454" ,"SRS3983999.pilon_04455" ,"SRS3983999.pilon_04456" ,"SRS3983999.pilon_04457" ,"SRS3983999.pilon_04458" ,"SRS3983999.pilon_04459" ,"SRS3983999.pilon_04460" ,"SRS3983999.pilon_04461" ,"SRS3983999.pilon_04462" ,"SRS3983999.pilon_04463" ,"SRS3983999.pilon_04464" ,"SRS3983999.pilon_04465" ,"SRS3983999.pilon_04466" ,"SRS3983999.pilon_04467" ,"SRS3983999.pilon_04468" ,"SRS3983999.pilon_04469" ,"SRS3983999.pilon_04470" ,"SRS3983999.pilon_04471" ,"SRS3983999.pilon_04472" ,"SRS3983999.pilon_04473" ,"SRS3983999.pilon_04474" ,"SRS3983999.pilon_04475" ,"SRS3983999.pilon_04476" ,"SRS3983999.pilon_04477" ,"SRS3983999.pilon_04478" ,"SRS3983999.pilon_04479" ,"SRS3983999.pilon_04480" ,"SRS3983999.pilon_04481" ,"SRS3983999.pilon_04482" ,"SRS3983999.pilon_04483" ,"SRS3983999.pilon_04484" ,"SRS3983999.pilon_04485" ,"SRS3983999.pilon_04486" ,"SRS3983999.pilon_04487" ,"SRS3983999.pilon_04488" ,"SRS3983999.pilon_04489" ,"SRS3983999.pilon_04490" ,"SRS3983999.pilon_04491" ,"SRS3983999.pilon_04492" ,"SRS3983999.pilon_04493" ,"SRS3983999.pilon_04494" ,"SRS3983999.pilon_04495" ,"SRS3983999.pilon_04496" ,"SRS3983999.pilon_04497" ,"SRS3983999.pilon_04498" ,"SRS3983999.pilon_04499" ,"SRS3983999.pilon_04500" ,"SRS3983999.pilon_04501" ,"SRS3983999.pilon_04502" ,"SRS3983999.pilon_04503" ,"SRS3983999.pilon_04504" ,"SRS3983999.pilon_04505" ,"SRS3983999.pilon_04506" ,"SRS3983999.pilon_04507" ,"SRS3983999.pilon_04508" ,"SRS3983999.pilon_04509" ,"SRS3983999.pilon_04510" ,"SRS3983999.pilon_04511" ,"SRS3983999.pilon_04512" ,"SRS3983999.pilon_04513" ,"SRS3983999.pilon_04514" ,"SRS3983999.pilon_04515" ,"SRS3983999.pilon_04516" ,"SRS3983999.pilon_04517" ,"SRS3983999.pilon_04518" ,"SRS3983999.pilon_04519" ,"SRS3983999.pilon_04520" ,"SRS3983999.pilon_04521" ,"SRS3983999.pilon_04522" ,"SRS3983999.pilon_04523" ,"SRS3983999.pilon_04524" ,"SRS3983999.pilon_04525" ,"SRS3983999.pilon_04526" ,"SRS3983999.pilon_04527" ,"SRS3983999.pilon_04528" ,"SRS3983999.pilon_04529" ,"SRS3983999.pilon_04530" ,"SRS3983999.pilon_04531" ,"SRS3983999.pilon_04532" ,"SRS3983999.pilon_04533" ,"SRS3983999.pilon_04534" ,"SRS3983999.pilon_04535" ,"SRS3983999.pilon_04536" ,"SRS3983999.pilon_04537" ,"SRS3983999.pilon_04538" ,"SRS3983999.pilon_04539" ,"SRS3983999.pilon_04540" ,"SRS3983999.pilon_04541" ,"SRS3983999.pilon_04542" ,"SRS3983999.pilon_04543" ,"SRS3983999.pilon_04544" ,"SRS3983999.pilon_04545" ,"SRS3983999.pilon_04546" ,"SRS3983999.pilon_04547" ,"SRS3983999.pilon_04548" ,"SRS3983999.pilon_04549" ,"SRS3983999.pilon_04550" ,"SRS3983999.pilon_04551" ,"SRS3983999.pilon_04552" ,"SRS3983999.pilon_04553" ,"SRS3983999.pilon_04554" ,"SRS3983999.pilon_04555" ,"SRS3983999.pilon_04556" ,"SRS3983999.pilon_04557" ,"SRS3983999.pilon_04558" ,"SRS3983999.pilon_04559" ,"SRS3983999.pilon_04560" ,"SRS3983999.pilon_04561" ,"SRS3983999.pilon_04562" ,"SRS3983999.pilon_04563" ,"SRS3983999.pilon_04564" ,"SRS3983999.pilon_04565" ,"SRS3983999.pilon_04566" ,"SRS3983999.pilon_04567" ,"SRS3983999.pilon_04568" ,"SRS3983999.pilon_04569" ,"SRS3983999.pilon_04570" ,"SRS3983999.pilon_04571" ,"SRS3983999.pilon_04572" ,"SRS3983999.pilon_04573" ,"SRS3983999.pilon_04574" ,"SRS3983999.pilon_04575" ,"SRS3983999.pilon_04576" ,"SRS3983999.pilon_04577" ,"SRS3983999.pilon_04578" ,"SRS3983999.pilon_04579" ,"SRS3983999.pilon_04580" ,"SRS3983999.pilon_04581" ,"SRS3983999.pilon_04582" ,"SRS3983999.pilon_04583" ,"SRS3983999.pilon_04584" ,"SRS3983999.pilon_04585" ,"SRS3983999.pilon_04586" ,"SRS3983999.pilon_04587" ,"SRS3983999.pilon_04588" ,"SRS3983999.pilon_04589" ,"SRS3983999.pilon_04590" ,"SRS3983999.pilon_04591" ,"SRS3983999.pilon_04592" ,"SRS3983999.pilon_04593" ,"SRS3983999.pilon_04594" ,"SRS3983999.pilon_04595" ,"SRS3983999.pilon_04596" ,"SRS3983999.pilon_04597" ,"SRS3983999.pilon_04598" ,"SRS3983999.pilon_04599" ,"SRS3983999.pilon_04600" ,"SRS3983999.pilon_04601" ,"SRS3983999.pilon_04602" ,"SRS3983999.pilon_04603" ,"SRS3983999.pilon_04604" ,"SRS3983999.pilon_04605" ,"SRS3983999.pilon_04606" ,"SRS3983999.pilon_04607" ,"SRS3983999.pilon_04608" ,"SRS3983999.pilon_04609" ,"SRS3983999.pilon_04610" ,"SRS3983999.pilon_04611" ,"SRS3983999.pilon_04612" ,"SRS3983999.pilon_04613" ,"SRS3983999.pilon_04614" ,"SRS3983999.pilon_04615" ,"SRS3983999.pilon_04616" ,"SRS3983999.pilon_04617" ,"SRS3983999.pilon_04618" ,"SRS3983999.pilon_04619" ,"SRS3983999.pilon_04620" ,"SRS3983999.pilon_04621" ,"SRS3983999.pilon_04622" ,"SRS3983999.pilon_04623" ,"SRS3983999.pilon_04624" ,"SRS3983999.pilon_04625" ,"SRS3983999.pilon_04626" ,"SRS3983999.pilon_04627" ,"SRS3983999.pilon_04628" ,"SRS3983999.pilon_04629" ,"SRS3983999.pilon_04630" ,"SRS3983999.pilon_04631" ,"SRS3983999.pilon_04632" ,"SRS3983999.pilon_04633" ,"SRS3983999.pilon_04634" ,"SRS3983999.pilon_04635" ,"SRS3983999.pilon_04636" ,"SRS3983999.pilon_04637" ,"SRS3983999.pilon_04638" ,"SRS3983999.pilon_04639" ,"SRS3983999.pilon_04640" ,"SRS3983999.pilon_04641" ,"SRS3983999.pilon_04642" ,"SRS3983999.pilon_04643" ,"SRS3983999.pilon_04644" ,"SRS3983999.pilon_04645" ,"SRS3983999.pilon_04646" ,"SRS3983999.pilon_04647" ,"SRS3983999.pilon_04648" ,"SRS3983999.pilon_04649"]
    print("Number of seqs {}".format(len(seqs_oi)))
    seqs_1 = faster_fasta_searching(infile, seqs_oi)

    z = run_blast(seqs_1,seqs_1)

    z = run_mmseqs(seqs_1, seqs_1)

    outfile = write_random_data()
    z = read_data_csv(outfile)
    Df = read_data_pandas(outfile)
    y = loop_data_pandas(Df)

    infile = "1250Genomes_HardCorePrototype/prokka/SRS3983999.pilon/SRS3983999.pilon.ffn"
    seqs_oi = ['SRS3983999.pilon_04684',
    'SRS3983999.pilon_04685',
    'SRS3983999.pilon_04686',
    'SRS3983999.pilon_04687',
    'SRS3983999.pilon_04688',
    'SRS3983999.pilon_04689',
    'SRS3983999.pilon_04690',
    'SRS3983999.pilon_04691',
    'SRS3983999.pilon_04692',
    'SRS3983999.pilon_04693',
    'SRS3983999.pilon_04694']

    seqs_1 = faster_fasta_searching(infile, seqs_oi)
    #print(seqs_1[:5])

    seqs_2 = read_seqs(infile, seqs_oi)
    #print(seqs_2[:5])
    assert [s.seq for s in seqs_1 if s.id == "SRS3983999.pilon_04690"][0] == [s.seq for s in seqs_2 if s.id == "SRS3983999.pilon_04690"][0]
    assert len(seqs_1) == len(seqs_2)
    assert sorted([x.id for x in seqs_1]) == sorted([x.id for x in seqs_2])