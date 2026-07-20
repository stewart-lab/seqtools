#!/usr/bin/env python3
# Run fastp QC/trimming on a directory of paired-end fastq files.
# Use this to inspect read quality before running the full pipeline.
#
# Usage:
#   conda activate star_env
#   python run_fastp.py -i /path/to/fastqs -o /path/to/output

import os
import json
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='Run fastp on a directory of paired-end fastq files')
parser.add_argument('-i', '--fastq_dir', type=str, required=True,
                    help='Directory containing paired-end .fastq/.fq.gz/.fastq.gz files')
parser.add_argument('-o', '--output_dir', type=str, required=True,
                    help='Output directory (must exist; a timestamped subdirectory will be created)')
args = parser.parse_args()


def main():
    fastq_dir = args.fastq_dir
    output_dir = args.output_dir

    _check_inputs(fastq_dir, output_dir)

    # create timestamped output dir
    current_datetime = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    run_dir = os.path.join(output_dir, current_datetime + '___fastp')
    os.makedirs(run_dir, exist_ok=True)
    print(f"Output directory: {run_dir}")

    _run_fastp(fastq_dir, run_dir)
    _write_summary(run_dir)

    print(f"\nDone. Summary written to {os.path.join(run_dir, '_filtering_summary.txt')}")


def _check_inputs(fastq_dir: str, output_dir: str) -> None:
    if not os.path.isdir(output_dir):
        raise Exception("Output directory does not exist")
    if not os.path.isdir(fastq_dir):
        raise Exception("Fastq directory does not exist")

    fastq_files = [f for f in os.listdir(fastq_dir)
                   if f.endswith('.fastq') or f.endswith('.fq.gz') or f.endswith('.fastq.gz')]
    if not fastq_files:
        raise Exception("No fastq files found in the fastq directory")

    r1_fastqs = [f for f in fastq_files if ('_R1_' in f) or ('1.fq.gz' in f) or ('_R1.' in f)]
    r2_fastqs = [f for f in fastq_files if ('_R2_' in f) or ('2.fq.gz' in f) or ('_R2.' in f)]
    if not r1_fastqs:
        raise Exception("No R1 fastq files found")
    if not r2_fastqs:
        raise Exception("No R2 fastq files found")
    if len(r1_fastqs) != len(r2_fastqs):
        raise Exception("Unequal number of R1 and R2 fastq files")
    for r1 in r1_fastqs:
        if '_R1_' in r1:
            r2 = r1.replace('_R1_', '_R2_')
        elif '1.fq.gz' in r1:
            r2 = r1.replace('1.fq.gz', '2.fq.gz')
        elif '_R1.' in r1:
            r2 = r1.replace('_R1.', '_R2.')
        else:
            if r2 not in r2_fastqs:
                raise Exception(f"No matching R2 file found for {r1}")

    print(f"Found {len(r1_fastqs)} sample(s) to process.")


def _run_fastp(fastq_dir: str, run_dir: str) -> None:
    r1_files = [f for f in os.listdir(fastq_dir)
                if ('_R1_' in f and f.endswith('.fastq'))
                or ('1.fq' in f and f.endswith('.gz'))
                or ('_R1_' in f and f.endswith('.fastq.gz'))
                or ('_R1.' in f and f.endswith('.fastq.gz'))]
    r1_files.sort()

    for fastq1 in r1_files:
        if '_R1_' in fastq1 and fastq1.endswith('.fastq.gz'):
            fastq2 = fastq1.replace('_R1_', '_R2_')
            html = fastq1.replace('.fastq.gz', '.html').replace('_R1_', '').replace('001.html', '.html')
        elif '_R1.' in fastq1 and fastq1.endswith('.fastq.gz'):
            fastq2 = fastq1.replace('_R1.', '_R2.')
            html = fastq1.replace('.fastq.gz', '.html').replace('_R1', '').replace('001.html', '.html')
        elif '_R1_' in fastq1 and fastq1.endswith('.fastq'):
            fastq2 = fastq1.replace('_R1_', '_R2_')
            html = fastq1.replace('.fastq', '.html').replace('_R1_', '').replace('001.html', '.html')
        elif '1.fq.gz' in fastq1:
            fastq2 = fastq1.replace('1.fq.gz', '2.fq.gz')
            html = fastq1.replace('1.fq.gz', '.html')

        fastq1_path = os.path.join(fastq_dir, fastq1)
        fastq2_path = os.path.join(fastq_dir, fastq2)
        fastq1_out = os.path.join(run_dir, fastq1)
        fastq2_out = os.path.join(run_dir, fastq2)
        html_path = os.path.join(run_dir, f"report_{html}")
        json_path = os.path.join(run_dir, f"report_{html.replace('.html', '.json')}")

        print(f"Running fastp on {fastq1} / {fastq2} ...")
        os.system(f'fastp --in1 {fastq1_path} --in2 {fastq2_path} '
                  f'--out1 {fastq1_out} --out2 {fastq2_out} '
                  f'-h {html_path} -j {json_path}')


def _write_summary(run_dir: str) -> None:
    fastp_reports = sorted([f for f in os.listdir(run_dir) if f.endswith('.json')])
    if not fastp_reports:
        print("No fastp JSON reports found — skipping summary.")
        return

    progress_char = '█'
    summary_path = os.path.join(run_dir, '_filtering_summary.txt')

    with open(summary_path, 'w') as out:
        for report_file in fastp_reports:
            with open(os.path.join(run_dir, report_file)) as f:
                data = json.load(f)

            total_before = float(data['summary']['before_filtering']['total_reads'])
            total_after  = float(data['summary']['after_filtering']['total_reads'])
            low_quality  = float(data['filtering_result']['low_quality_reads'])
            pct_kept     = round((total_after / total_before) * 100, 2)

            sample_name = report_file.replace('.json', '')

            before_bar = progress_char * int(total_before / 2_000_000)
            after_bar  = progress_char * int(total_after  / 2_000_000)
            before_bar = ' '.join(before_bar[i:i+10] for i in range(0, len(before_bar), 10))
            after_bar  = ' '.join(after_bar[i:i+10]  for i in range(0, len(after_bar),  10))

            out.write(f"Sample: {sample_name}\n")
            out.write(f"Total reads before filtering: {total_before}\n")
            out.write(f"Total reads after filtering:  {total_after}\n")
            out.write(f"Low quality reads:            {low_quality}\n")
            out.write(f"Percent of reads kept:        {pct_kept}%\n")
            out.write(f"Bar chart before filtering:   ---> {before_bar}\n")
            out.write(f"Bar chart after filtering:    ---> {after_bar}\n")
            out.write("\n\n")

            # also echo to stdout
            print(f"  {sample_name}: {pct_kept}% reads kept ({total_after:.0f} / {total_before:.0f})")


if __name__ == '__main__':
    main()
