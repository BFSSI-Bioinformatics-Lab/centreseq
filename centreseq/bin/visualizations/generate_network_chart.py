from pathlib import Path

NETWORK_GRAPH_HTML_TEMPLATE = Path(__file__).parents[1] / 'visualizations' / 'network_graph_template.html'
STATIC_DIR = Path(__file__).parents[1] / 'visualizations' / 'static'
assert NETWORK_GRAPH_HTML_TEMPLATE.exists()
assert STATIC_DIR.exists()


def generate_network_chart(pairwise_gene_count_report: Path, roary_report: Path, network_coding: Path, outdir: Path):
    # Read in the file
    with open(NETWORK_GRAPH_HTML_TEMPLATE, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('ROARY_GENE_COUNT_REPORT_VARIABLE', f"/reports/{roary_report.name}")
    filedata = filedata.replace('NETWORK_CODING', network_coding.name)
    filedata = filedata.replace('PAIRWISE_GENE_COUNT_REPORT_TSV', f"reports/{pairwise_gene_count_report.name}")
    filedata = filedata.replace('STATIC_DIR', str(STATIC_DIR))

    # Write the file out again
    outfile = outdir / 'network_graph.html'
    if outfile.exists():
        outfile.unlink()
    with open(str(outfile), 'w') as file:
        file.write(filedata)

    # Symlink the static files to static directory in outdir
    staticdir = outdir / 'static'
    staticdir.mkdir(exist_ok=True)
    staticfiles = list(STATIC_DIR.glob("*"))
    for s in staticfiles:
        s_ = staticdir / s.name
        if s_.exists():
            continue
        s_.symlink_to(s)

    return outfile


def generate_network_chart_coding_file(outdir: Path, sample_id_list: list) -> Path:
    network_coding = outdir / 'network_graph_coding.tsv'
    if network_coding.exists():
        network_coding.unlink()
    with open(str(network_coding), 'w') as file:
        file.write(f"sample_id\tgroup_id\n")
        for sample in sample_id_list:
            file.write(f"{sample}\tdefault\n")
    return network_coding
