import os
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='RNA-Seq pipeline')
parser.add_argument('-i', '--fastq_dir', type=str, help='Directory containing input fastq files', required=True)
parser.add_argument('-r', '--rsem_reference_dir', type=str, help='Path to RSEM reference', required=True)
parser.add_argument('-o', '--output_dir', type=str, help='Output directory', required=True)
args = parser.parse_args()

def main():
    fastq_dir = args.fastq_dir
    output_dir = args.output_dir
    rsem_reference = args.rsem_reference_dir
    output_dir = _make_output_dir(output_dir)

    # TODO: do a dry run for the input fastqs to make sure everything is good?
    # in general probably should have more error checking, i.e. make sure inputs are paired end, etc.

    # TODO: do a 'conda list' and write to a README file in the output directory
    # need to include flags

    # TODO: output genome version as well

    # run fastqc on input fastqs
    # fastqc1_dir = os.path.join(output_dir, '001_fastqc')
    # _run_fastqc(fastq_dir, fastqc1_dir)

    # trim input fastqs
    fastp_dir = os.path.join(output_dir, '002_fastq_trimmed')
    _trim_fastqs(fastq_dir, fastp_dir)

    # run fastqc on trimmed fastqs
    # fastqc2_dir = os.path.join(output_dir, '003_fastqc')
    # _run_fastqc(fastp_dir, fastqc2_dir)

    # run rsem on trimmed fastqs
    rsem_dir = os.path.join(output_dir, '004_rsem')
    _run_rsem(fastp_dir, rsem_dir, rsem_reference)

    # get counts matrix from RSEM output
    counts_dir = os.path.join(output_dir, '005_count_matrix')
    _get_counts_csv(rsem_dir, counts_dir)

    # TODO: make contrasts and such for DESeq2?


def _run_fastqc(input_dir: str, output_dir: str) -> None:
    os.makedirs(output_dir)
    fastqs = [f for f in os.listdir(input_dir) if f.endswith('.fastq')]
    for fastq in fastqs:
        # get input paths
        fastq_path = os.path.join(input_dir, fastq)

        # run fastqc
        os.system(f'fastqc --outdir={output_dir} {fastq_path}')

def _trim_fastqs(input_dir: str, output_dir: str) -> None:
    os.makedirs(output_dir)
    fastqs = [f for f in os.listdir(input_dir) if '_R1_' in f and f.endswith('.fastq')]
    for fastq1 in fastqs:
        # get input paths
        fastq2 = fastq1.replace('_R1_', '_R2_')
        html = fastq1.replace('.fastq', '.html').replace('_R1_', '').replace('001.html', '.html')
        fastq1_path = os.path.join(input_dir, fastq1)
        fastq2_path = os.path.join(input_dir, fastq2)

        # get output paths
        fastq1_output_path = os.path.join(output_dir, fastq1)
        fastq2_output_path = os.path.join(output_dir, fastq2)
        html_path = os.path.join(output_dir, html)

        # run fastp
        os.system(f'fastp --in1 {fastq1_path} --in2 {fastq2_path} --out1 {fastq1_output_path} --out2 {fastq2_output_path} -h {html_path}')
    
    # zip up the .html files for download
    # os.system(f'zip -r {output_dir}/__fastp_html_reports.zip {output_dir}/*.html')

def _run_rsem(input_dir: str, output_dir: str, rsem_reference: str) -> None:
    os.makedirs(output_dir)

    # change working directory... fixes super annoying rsem bug.
    # RSEM tries to create a folder called {sample_name}.stat in the current working directory.
    # but in a docker container the cwd is not writable by default.
    # so we change the cwd to the output directory, which is presumably writable.
    # https://github.com/deweylab/RSEM/blob/8bc1e2115493c0cdf3c6bee80ef7a21a91b2acce/rsem-calculate-expression#L372
    original_wd = os.getcwd()
    os.chdir(output_dir)

    fastqs = [f for f in os.listdir(input_dir) if '_R1_' in f and f.endswith('.fastq')]
    rsem_temporary_dir = os.path.join(output_dir, 'rsem_temp')
    for fastq1 in fastqs:
        fastq2 = fastq1.replace('_R1_', '_R2_')
        sample_name = fastq1.split('_R1_')[0]
        fastq1_path = os.path.join(input_dir, fastq1)
        fastq2_path = os.path.join(input_dir, fastq2)
        os.system(f'rsem-calculate-expression --num-threads 40 --star --temporary-folder {rsem_temporary_dir} --paired-end {fastq1_path} {fastq2_path} {rsem_reference} {sample_name}')

    # change working directory back to original working directory.
    # probably not necessary, but still.
    os.chdir(original_wd)

def _get_counts_csv(rsem_dir: str, output_dir: str) -> str:
    os.makedirs(output_dir)
    # build the counts.csv file from the output in the RSEM directory
    os.system(f'Rscript /src/tximport.R --rsem_dir {rsem_dir} --output_dir {output_dir}')

def _make_output_dir(output_dir: str) -> str:
    current_datetime = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    output_dir = os.path.join(output_dir, current_datetime + "___rnatools")
    os.makedirs(output_dir)
    return output_dir

if __name__ == '__main__':
    main()