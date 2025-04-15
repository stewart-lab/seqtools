import os
import argparse
import math
import csv
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--counts', type=str, help='Path to counts.csv', required=True)
parser.add_argument('--design', type=str, help='Path to design.csv', required=True)
parser.add_argument('--contrasts', type=str, help='Path to contrasts.csv', required=True)
parser.add_argument('--padj', type=float, help='Filter DESeq2 results by padj', default=0.05)
parser.add_argument('--lfc', type=float, help='Filter DESeq2 results by log2 fold change', default=0.5)
parser.add_argument('--gmt', type=str, help='Path to a gene set .gmt file to use for pathway analysis', required=True)
parser.add_argument('--extra_plots', type=str, help='Path to a csv file with extra plots to make', default=None)
parser.add_argument('-o', '--output', type=str, help='Path to output directory. Directory must exist.', required=True)
args = parser.parse_args()

def main():
    # parse args
    counts_path = args.counts
    design_path = args.design
    contrasts_path = args.contrasts
    output_dir = args.output
    max_padj = args.padj
    min_lfc = args.lfc
    gmt_path = args.gmt
    extra_plots_csv = args.extra_plots
    min_max = 10  # the minimum max value to keep the gene (filters out lowly expressed genes)

    # check for valid inputs
    _check_inputs(counts_path, design_path, contrasts_path, output_dir)

    # construct output directory paths
    deseq2_dir = os.path.join(output_dir, '005_deseq_output')
    filtered_deseq_dir = os.path.join(output_dir, '006_deseq_output_filtered')
    sig_genes_dir = os.path.join(output_dir, '007_differentially_expressed_genes')
    pathway_analysis_dir = os.path.join(output_dir, '008_pathway_analysis')
    plots_dir = os.path.join(output_dir, '009_plots')

    # clear anything in the output folders so we don't accidentally look at old results
    os.system(f'rm -rf {deseq2_dir}/*')
    os.system(f'rm -rf {filtered_deseq_dir}/*')
    os.system(f'rm -rf {sig_genes_dir}/*')
    os.system(f'rm -rf {pathway_analysis_dir}/*')
    os.system(f'rm -rf {plots_dir}/*')

    # run the pipeline
    _run_deseq(counts_path, design_path, contrasts_path, deseq2_dir)
    _filter_deseq_results(deseq2_dir, filtered_deseq_dir, sig_genes_dir, min_max, max_padj, min_lfc)
    _run_gprofiler(sig_genes_dir, pathway_analysis_dir, gmt_path)
    _make_qc_plots(filtered_deseq_dir, plots_dir)
    _make_extra_plots(extra_plots_csv, filtered_deseq_dir, gmt_path, plots_dir, max_padj, min_lfc)
    _make_contrast_plots(filtered_deseq_dir, pathway_analysis_dir, gmt_path, plots_dir, max_padj, min_lfc)


### pipeline steps
def _check_inputs(counts_path: str, design_path: str, contrasts_path: str, output_dir: str):
    print("checking inputs...")

    # check that counts, design, and contrasts exist and they are .csv files
    if not os.path.exists(counts_path):
        raise FileNotFoundError(f'Counts file not found at {counts_path}')
    
    if not os.path.exists(design_path):
        raise FileNotFoundError(f'Design file not found at {design_path}')
    
    if not os.path.exists(contrasts_path):
        raise FileNotFoundError(f'Contrasts file not found at {contrasts_path}')
    
    if not counts_path.endswith('.csv'):
        raise ValueError(f'Counts file must be .csv')
    
    if not design_path.endswith('.csv'):
        raise ValueError(f'Design file must be .csv')
    
    if not contrasts_path.endswith('.csv'):
        raise ValueError(f'Contrasts file must be .csv')
    
    # check that output directory exists
    if not os.path.exists(output_dir):
        raise FileNotFoundError(f'Output directory not found at {output_dir}')
    
    print('Inputs are valid!')

def _run_deseq(counts_path: str, design_path: str, contrasts_path: str, deseq2_dir: str) -> None:
    print("running DESeq2...")

    # make output dir
    os.makedirs(deseq2_dir, exist_ok=True)

    # run deseq2
    this_dir_path = os.path.dirname(os.path.realpath(__file__))
    os.system(f'Rscript {this_dir_path}/deseq2.r --counts {counts_path} --design {design_path} --contrasts {contrasts_path} --output_dir {deseq2_dir}')

def _filter_deseq_results(deseq2_dir: str, filtered_deseq_dir: str, sig_genes_dir, min_max: float, max_padj: float, min_lfc: float) -> None:
    print("filtering DESeq2 results...")

    # make output dir
    os.makedirs(filtered_deseq_dir, exist_ok=True)

    # filter contrasts
    deseq_output_files = [f for f in os.listdir(deseq2_dir) if f.endswith('.csv')]
    for file in deseq_output_files:
        # skip anything that isn't a contrast (e.g. normalized_counts.csv)
        is_contrast = '_vs_' in file
        file_path = os.path.join(deseq2_dir, file)
        if not is_contrast:
            continue

        # read deseq contrast result
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # write filtered deseq results to new file (_filtered.csv)
        filtered = _filter_deseq_lines(lines, 1, 0, min_max, 3, 1)
        filtered_file_path = os.path.join(filtered_deseq_dir, file.replace('.csv', '_filtered.csv'))
        with open(filtered_file_path, 'w') as f:
            f.writelines(filtered)

        # write significant deseq results to new file (_significant.csv)
        de = _filter_deseq_lines(lines, max_padj, min_lfc, min_max, 3, 1)
        de_file_path = os.path.join(sig_genes_dir, file.replace('.csv', '_significant.csv'))
        with open(de_file_path, 'w') as f:
            f.writelines(de)

    # filter normalized_counts.csv
    normalized_counts_path = os.path.join(deseq2_dir, 'normalized_counts.csv')
    if os.path.exists(normalized_counts_path):
        with open(normalized_counts_path, 'r') as f:
            lines = f.readlines()

        filtered = _filter_deseq_lines(lines, 1, 0, min_max, 3, 1)

        # write filtered deseq2 results to new file (_filtered.csv)
        filtered_file_path = os.path.join(filtered_deseq_dir, 'normalized_counts_filtered.csv')
        with open(filtered_file_path, 'w') as f:
            f.writelines(filtered)

def _run_gprofiler(significant_deseq_output_dir: str, pathway_analysis_dir: str, gmt_path: str) -> None:
    print("running pathway analysis...")
    this_dir_path = os.path.dirname(os.path.realpath(__file__))
    n_required_de_genes = 50
    max_term_size = 300

    # upload GMT and get the custom organism token
    result = subprocess.run(
        ["Rscript", f"{this_dir_path}/gprofiler_upload_gmt.R", "--gmt", gmt_path],
        capture_output=True,        # Capture standard output and error
        text=True                   # Decode the output as a string
    )
    organism = result.stdout.strip()

    if not organism:
        raise ValueError(f'Error uploading custom .gmt file to gProfiler. The server may not be available - try the web version.')

    print(f'**** Uploaded custom .gmt to gProfiler. Organism token: {organism}')

    # make output dir
    os.makedirs(pathway_analysis_dir, exist_ok=True)

    # run gprofiler on each filtered deseq2 result
    deseq_results = [f for f in os.listdir(significant_deseq_output_dir) if f.endswith('.csv')]
    for deseq_result in deseq_results:
        output_filename = deseq_result.replace('_significant.csv', '_pathways.csv')
        output_filepath = os.path.join(pathway_analysis_dir, output_filename)
        deseq_result_path = os.path.join(significant_deseq_output_dir, deseq_result)

        with open(deseq_result_path, 'r') as f:
            lines = f.readlines()

        # if there are not enough DE genes to run gProfiler, skip this contrast
        n_de_genes = len(lines) - 1
        if n_de_genes < n_required_de_genes:
            with open(output_filepath, 'w') as f:
                f.write(f'Not enough differentially expressed genes to run pathway analysis! Needed at least {n_required_de_genes} but only had {n_de_genes}')
            continue

        # run gprofiler
        os.system(f'Rscript {this_dir_path}/gprofiler.r --genes {deseq_result_path} --output_dir {pathway_analysis_dir} --organism {organism} --mthreshold {max_term_size}')

def _make_qc_plots(filtered_deseq_dir: str, plots_dir: str) -> None:
    # make output dir
    os.makedirs(plots_dir, exist_ok=True)

    # read normalized count matrix
    normalized_counts_path = os.path.join(filtered_deseq_dir, 'normalized_counts_filtered.csv')
    if os.path.exists(normalized_counts_path):
        with open(normalized_counts_path, 'r') as f:
            lines = f.readlines()
    else:
        print('No normalized_counts.csv file found. Skipping QC plots.')
        return
    
    ######### experiment-level plots #########
    # PCA
    # clustermap (+sigclust2?)

def _make_contrast_plots(filtered_deseq_dir: str, pathway_analysis_dir: str, gmt_path: str, plots_dir: str, max_padj: float, min_lfc: float) -> None:
    # get the DESeq2 contrast results
    deseq_results = [f for f in os.listdir(filtered_deseq_dir) if f.endswith('.csv') and '_vs_' in f]
    known_genesets = _read_geneset_gmt(gmt_path)
    for deseq_result in deseq_results:
        contrast_name = deseq_result.replace('_filtered.csv', '')
        deseq_result_path = os.path.join(filtered_deseq_dir, deseq_result)

        with open(deseq_result_path, 'r') as f:
            deseq_lines = f.readlines()

        # make a dir for this contrast
        contrast_dir = os.path.join(plots_dir, contrast_name)
        os.makedirs(contrast_dir, exist_ok=True)

        ######### contrast-level plots #########
        # make volcano plot of all genes for this contrast
        output_file_path = os.path.join(contrast_dir, f'{contrast_name}_VOLCANO')
        _plot_volcano(contrast_name, deseq_result_path, output_file_path, max_padj, min_lfc)

        # GSVA plot
        # TODO: should this be on significant genes only?
        # if os.path.exists(gmt_path) and gmt_path.endswith('.gmt'):
        #     gsva_file_path = os.path.join(contrast_dir, f'{contrast_name}_GSVA')
        #     _plot_gsva(deseq_result_path, gmt_path, gsva_file_path)

        # TODO: make MA plot
        # TODO: make clustermap


        ######### geneset-level plots #########
        contrast_pathways_file = os.path.join(pathway_analysis_dir, f'{contrast_name}_pathways.csv')
        significant_pathways = _get_gene_sets_from_pathway_analysis(contrast_pathways_file, known_genesets)

        for geneset_name, geneset_members in significant_pathways.items():
            deseq_geneset_lines = _filter_deseq_lines(deseq_lines, geneset=geneset_members)

            sanitized_geneset_name = ''.join(e for e in geneset_name if e.isalnum() or e in [' ', '_', '-']).replace(' ', '_')
            geneset_dir = os.path.join(contrast_dir, "pathways", sanitized_geneset_name)
            os.makedirs(geneset_dir, exist_ok=True)

            csv_file_path = os.path.join(geneset_dir, f'{sanitized_geneset_name}.csv')
            heatmap_file_path = os.path.join(geneset_dir, f'{sanitized_geneset_name}_HEATMAP')
            volcano_file_path = os.path.join(geneset_dir, f'{sanitized_geneset_name}_VOLCANO')

            # write DESeq results for this gene set to .csv
            with open(csv_file_path, 'w') as f:
                for line in deseq_geneset_lines:
                    f.write(line)

            plot_title = f'{contrast_name} -- {geneset_name}'.replace('_', ' ')
            _plot_heatmap(plot_title, csv_file_path, heatmap_file_path)
            _plot_volcano(plot_title, csv_file_path, volcano_file_path, max_padj, min_lfc)

def _make_extra_plots(extra_plots_csv: str, filtered_deseq_dir: str, gmt_path: str, plots_dir: str, max_padj: float, min_lfc: float) -> None:
    valid_plot_types = ['heatmap', 'volcano', 'gsva']
    valid_gene_prefixes = ['genes', 'pathways']
    valid_samples_prefixes = ['samples', 'contrast']
    
    with open(extra_plots_csv, 'r') as f:
        extra_plot_lines = f.readlines()

    header = extra_plot_lines[0]
    spl_header = [x.strip() for x in header.split(',')]
    plot_type_idx = spl_header.index('plot')
    genes_idx = spl_header.index('genes')
    samples_idx = spl_header.index('samples')
    plot_title_idx = spl_header.index('title')

    for extra_plot_line in extra_plot_lines[1:]:
        spl_extra_plot_line = extra_plot_line.strip().split(',')
        plot_type = spl_extra_plot_line[plot_type_idx]
        genes_str = spl_extra_plot_line[genes_idx]
        samples_str = spl_extra_plot_line[samples_idx]
        plot_title = spl_extra_plot_line[plot_title_idx]

        if genes_str.count(':') != 1:
            print(f'WARNING: Invalid genes string "{genes_str}" in {extra_plots_csv}. Must be in the format "genes:gene1;gene2;gene3" or "pathways:gene_set_name".')
            continue

        if samples_str.count(':') != 1:
            print(f'WARNING: Invalid samples string "{samples_str}" in {extra_plots_csv}. Must be in the format "samples:sample1;sample2;sample3" or "contrast:contrast_name".')
            continue

        if plot_type not in valid_plot_types:
            print(f'WARNING: Invalid plot type "{plot_type}" in {extra_plots_csv}. Valid types are {valid_plot_types}.')
            continue

        genes_prefix = genes_str.split(':')[0]
        samples_prefix = samples_str.split(':')[0]

        if genes_prefix not in valid_gene_prefixes:
            print(f'WARNING: Invalid genes prefix "{genes_prefix}" in {extra_plots_csv}. Valid prefixes are {valid_gene_prefixes}.')
            continue

        if samples_prefix not in valid_samples_prefixes:
            print(f'WARNING: Invalid samples prefix "{samples_prefix}" in {extra_plots_csv}. Valid prefixes are {valid_samples_prefixes}.')
            continue

        if genes_prefix == 'genes':
            genes = genes_str.split(':')[1].split(';')
        elif genes_prefix == 'pathways':
            pathway_names = genes_str.split(':')[1].split(';')
            known_genesets = _read_geneset_gmt(gmt_path)
            genes = []
            for pathway_name in pathway_names:
                if pathway_name not in known_genesets:
                    print(f'WARNING: Gene set "{pathway_name}" not found in known gene sets. Skipping extra plot.')
                    continue
                genes += known_genesets[pathway_name]

        if samples_prefix == 'samples':
            samples = samples_str.split(':')[1].split(';')

            # read the normalized_counts_filtered.csv file. filter by genes and sample names
            normalized_counts_path = os.path.join(filtered_deseq_dir, 'normalized_counts_filtered.csv')
            if not os.path.exists(normalized_counts_path):
                print(f'WARNING: normalized_counts_filtered.csv not found at {normalized_counts_path}. Skipping extra plot.')
                continue
            with open(normalized_counts_path, 'r') as f:
                lines = f.readlines()
                # filter by samples
            spl_header = lines[0].strip().split(',')
            sample_idxs = [spl_header.index(x.strip()) for x in samples]
            lines = [','.join([line.split(',')[0].strip()] + [line.split(',')[i].strip() for i in sample_idxs]) for line in lines]
            lines = [x + '\n' for x in lines]
        elif samples_prefix == 'contrast':
            contrast_name = samples_str.split(':')[1]
            deseq_result_path = os.path.join(filtered_deseq_dir, f'{contrast_name}_filtered.csv')
            if not os.path.exists(deseq_result_path):
                print(f'WARNING: Contrast file not found at {deseq_result_path}. Skipping extra plot.')
                continue
            with open(deseq_result_path, 'r') as f:
                deseq_result_header = f.readline()
            spl_deseq_result_header = deseq_result_header.split(',')
            padj_idx = spl_deseq_result_header.index('padj')
            samples = spl_deseq_result_header[padj_idx + 1:]
            with open(deseq_result_path, 'r') as f:
                lines = f.readlines()

        # filter by gene (if not doing GSVA)
        if plot_type != 'gsva':
            filtered_lines = _filter_deseq_lines(lines, geneset=genes)
        else:
            filtered_lines = lines

        

        # determine what folder to write the plot+data to
        sanitized_plot_title = ''.join(e for e in plot_title if e.isalnum() or e in [' ', '_', '-']).replace(' ', '_')
        write_dir = plots_dir
        if samples_prefix == 'contrast':
            write_dir = os.path.join(write_dir, contrast_name, 'extra', sanitized_plot_title)
            os.makedirs(write_dir, exist_ok=True)
        else:
            write_dir = os.path.join(plots_dir, 'extra', sanitized_plot_title)
            os.makedirs(write_dir, exist_ok=True)

        # write filtered deseq results to new file
        data_csv_path = os.path.join(write_dir, f'{sanitized_plot_title}.csv')
        with open(data_csv_path, 'w') as f:
            f.writelines(filtered_lines)

        # make the plot
        if plot_type == 'heatmap':
            plot_path = os.path.join(write_dir, f'{sanitized_plot_title}_HEATMAP')
            _plot_heatmap(plot_title, data_csv_path, plot_path)
        elif plot_type == 'volcano':
            plot_path = os.path.join(write_dir, f'{sanitized_plot_title}_VOLCANO')
            _plot_volcano(plot_title, data_csv_path, plot_path, max_padj, min_lfc)
        elif plot_type == 'gsva':
            plot_path = os.path.join(write_dir, f'{sanitized_plot_title}_GSVA')
            _plot_gsva(data_csv_path, gmt_path, plot_path, pathway_names)
            

### helpers
def _filter_deseq_lines(deseq_lines: list, max_padj: float = math.inf, min_lfc: float = 0, min_max: float = 0, min_max_count: int = 1, max_percent_zero: float = 1, geneset: list = []) -> list:
    header = deseq_lines[0]
    spl_header = [x.strip() for x in header.split(',')]
    gene_idx = spl_header.index('gene')

    if 'padj' in spl_header:
        # contrast
        padj_idx = spl_header.index('padj')
        lfc_idx = spl_header.index('log2FoldChange')
        mmt_idxs = [spl_header.index(x) for x in spl_header[padj_idx + 1:]]
    else:
        # normalized count matrix
        padj_idx = -1
        lfc_idx = -1
        mmt_idxs = [spl_header.index(x) for x in spl_header[1:]]

    filtered = [header]

    for line in deseq_lines[1:]:
        spl_line = line.split(',')
        
        # filter by gene name
        gene = spl_line[gene_idx]
        if geneset and gene not in geneset:
            continue

        # filter by padj. keeps NAs if min_padj is inf
        if padj_idx >= 0 and not math.isinf(max_padj):
            padj_str = spl_line[padj_idx]
            if padj_str == 'NA':
                continue
            if float(padj_str) > max_padj:
                continue

        # filter by lfc
        if lfc_idx >= 0:
            lfc = float(spl_line[lfc_idx])
            if abs(lfc) < min_lfc:
                continue

        mmts = [float(spl_line[idx]) for idx in mmt_idxs]

        # filter by max value
        mmts_greater_than_minmax = [x for x in mmts if x >= min_max]
        if len(mmts_greater_than_minmax) < min_max_count:
            continue

        # filter by percent zero
        percent_zero = mmts.count(0) / len(mmts)
        if percent_zero > max_percent_zero:
            continue

        filtered.append(line)

    if padj_idx >= 0:
        # sort by padj, ascending. put NAs on the bottom
        filtered = [filtered[0]] + sorted(filtered[1:], key=lambda x: float(x.split(',')[padj_idx]) if x.split(',')[padj_idx] != 'NA' else 1.1)

    return filtered

def _get_gene_sets_from_pathway_analysis(pathway_analysis_file_path: str, known_genesets: dict) -> dict:
    gene_sets = dict()

    if not os.path.exists(pathway_analysis_file_path):
        print(f'WARNING: Pathway analysis file not found at {pathway_analysis_file_path}.')
        return gene_sets

    with open(pathway_analysis_file_path, 'r') as f:
        csvreader = csv.reader(f)
        spl_lines = [x for x in csvreader]

    if len(spl_lines) <= 1:
        return gene_sets

    spl_gprofiler_header = spl_lines[0]
    term_name_idx = spl_gprofiler_header.index('term_name')

    for spl_line in spl_lines[1:]:
        term_name = spl_line[term_name_idx].strip()
        term_members = known_genesets.get(term_name, [])
        if not term_members:
            print(f'WARNING: Gene set "{term_name}" not found in known gene sets. Skipping.')
            continue
        gene_sets[term_name] = term_members

    return gene_sets    

def _read_geneset_gmt(gmt_file_path: str) -> dict:
    genesets = dict()

    if not os.path.exists(gmt_file_path):
        print(f'WARNING: GMT file not found at {gmt_file_path}.')
        return genesets
    
    with open(gmt_file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        spl_line = line.split('\t')
        geneset_name = spl_line[0]
        geneset_description = spl_line[1]
        geneset_genes = [x.strip() for x in spl_line[2:] if x.strip()]
        genesets[geneset_name] = geneset_genes

    return genesets

### plotting
def _plot_volcano(plot_title: str, deseq_result_path: str, output_file_path: str, min_padj: float, min_lfc: float) -> None:
    # check if csv exists
    if not os.path.exists(deseq_result_path):
        print(f'WARNING: DESeq2 result file not found at {deseq_result_path}. Skipping volcano plot.')
        return

    this_dir_path = os.path.dirname(os.path.realpath(__file__))
    os.system(f'Rscript {this_dir_path}/plotVolcano.R --data {deseq_result_path} --output {output_file_path} --title "{plot_title}" --padj_cutoff {min_padj} --l2fc_cutoff {min_lfc}')

def _plot_heatmap(plot_title: str, deseq_result_path: str, output_file_path: str) -> None:
    # check if csv exists
    if not os.path.exists(deseq_result_path):
        print(f'WARNING: CSV file not found at {deseq_result_path}. Skipping heatmap plot.')
        return

    this_dir_path = os.path.dirname(os.path.realpath(__file__))
    os.system(f'Rscript {this_dir_path}/plotHeatmap.R --data {deseq_result_path} --output {output_file_path} --title "{plot_title}" --cluster_rows TRUE --log2 TRUE --hide_row_names FALSE --subtract_row_means TRUE --legend_title "Relative Expression"')

def _plot_gsva(deseq_result_path: str, gmt_path: str, output_file_path: str, geneset_filter: list = []) -> None:
    # check if csv exists
    if not os.path.exists(deseq_result_path):
        print(f'WARNING: CSV file not found at {deseq_result_path}. Skipping heatmap plot.')
        return
    if not os.path.exists(gmt_path):
        print(f'WARNING: GMT file not found at {gmt_path}. Skipping heatmap plot.')
        return
    
    # TEMP: if output_file_path does not end in .png, add it
    if not output_file_path.endswith('.png'):
        output_file_path += '.png'

    this_dir_path = os.path.dirname(os.path.realpath(__file__))

    if geneset_filter:
        os.system(f'Rscript {this_dir_path}/gsva.R --data {deseq_result_path} --gene_set {gmt_path} --output {output_file_path} --filter {",".join(geneset_filter)}')
    else:
        os.system(f'Rscript {this_dir_path}/gsva.R --data {deseq_result_path} --gene_set {gmt_path} --output {output_file_path}')

if __name__ == '__main__':
    main()