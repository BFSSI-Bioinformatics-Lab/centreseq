from pathlib import Path

NETWORK_GRAPH_HTML_TEMPLATE = Path(__file__).parents[1] / 'visualizations' / 'network_graph_template.html'
# D3_V5 = Path(__file__).parents[1] / 'visualizations' / 'static' / 'd3.v5.min.js'
D3_V5 = 'https://d3js.org/d3.v5.min.js'
# D3_SIMPLE_SLIDER = Path(__file__).parents[1] / 'visualizations' / 'static' / 'd3-simple-slider.js'
D3_SIMPLE_SLIDER = 'https://unpkg.com/d3-simple-slider'
# print(NETWORK_GRAPH_HTML_TEMPLATE)
assert NETWORK_GRAPH_HTML_TEMPLATE.exists()


# D3_V5 = '/graph_assets/d3.v5.min.js'
# D3_SIMPLE_SLIDER = './graph_assets/d3-simple-slider.js'

def generate_network_chart(pairwise_gene_count_report: Path, roary_report: Path, outdir: Path):
    # Read in the file
    with open(NETWORK_GRAPH_HTML_TEMPLATE, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('ROARY_GENE_COUNT_REPORT_VARIABLE', roary_report.name)
    filedata = filedata.replace('PAIRWISE_GENE_COUNT_REPORT_TSV', pairwise_gene_count_report.name)
    filedata = filedata.replace('D3_V5', str(D3_V5))
    filedata = filedata.replace('D3_SIMPLE_SLIDER', str(D3_SIMPLE_SLIDER))

    # Write the file out again
    outfile = outdir / 'network_graph.html'
    if outfile.exists():
        outfile.unlink()
    with open(str(outfile), 'w') as file:
        file.write(filedata)

    return outfile
