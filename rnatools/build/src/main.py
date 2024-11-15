import os
import argparse
import json
from datetime import datetime

parser = argparse.ArgumentParser(description='RNA-Seq pipeline')
parser.add_argument('-i', '--fastq_dir', type=str, help='Directory containing input fastq files', required=True)
parser.add_argument('-r', '--reference_dir', type=str, help='Directory containing the .gtf and .fa genome reference files', required=True)
parser.add_argument('-o', '--output_dir', type=str, help='Output directory', required=True)
args = parser.parse_args()

def main():
    fastq_dir = args.fastq_dir
    output_dir = args.output_dir
    reference_dir = args.reference_dir

    # TODO: do a dry run for the input fastqs to make sure everything is good?
    _check_inputs(fastq_dir, reference_dir, output_dir)

    output_dir = _make_output_dir(output_dir)

    # TODO: do a 'conda list' and write to a README file in the output directory
    # need to include flags
    # TODO: output genome version as well

    # trim input fastqs
    fastp_dir = os.path.join(output_dir, '001_fastq_trimmed')
    _trim_fastqs(fastq_dir, fastp_dir)

    # summarize fastp reports
    _summarize_fastp_reports(fastp_dir, output_dir)

    # run rsem on trimmed fastqs
    rsem_dir = os.path.join(output_dir, '002_rsem')
    _run_rsem(fastp_dir, rsem_dir, reference_dir)

    # summarize fastp reports again, this time with RSEM counts
    _summarize_fastp_reports(fastp_dir, output_dir, rsem_dir)

    # get counts matrix from RSEM output
    counts_dir = os.path.join(output_dir, '003_count_matrix')
    _write_counts_csv(rsem_dir, counts_dir)

    # TODO: make contrasts and such for DESeq2?

def _check_inputs(fastq_dir: str, reference_dir: str, output_dir: str) -> None:
    if not os.path.isdir(fastq_dir):
        raise Exception("Fastq directory does not exist")
    
    if not os.path.isdir(reference_dir):
        raise Exception("Reference directory does not exist")
    
    # check to see if the fastq dir contains any fastq files
    fastq_files = [f for f in os.listdir(fastq_dir) if f.endswith('.fastq')]
    if not fastq_files:
        raise Exception("No .fastq files found in the fastq directory")
    
    # all fastq files should be paired end
    r1_fastqs = [f for f in fastq_files if '_R1_' in f]
    r2_fastqs = [f for f in fastq_files if '_R2_' in f]
    if not r1_fastqs:
        raise Exception("No _R1_ fastq files found in the fastq directory")
    if not r2_fastqs:
        raise Exception("No _R2_ fastq files found in the fastq directory")
    if len(r1_fastqs) != len(r2_fastqs):
        raise Exception("Unequal number of _R1_ and _R2_ fastq files found in the fastq directory")
    for r1 in r1_fastqs:
        r2 = r1.replace('_R1_', '_R2_')
        if r2 not in r2_fastqs:
            raise Exception(f"Paired end _R2_ fastq file not found for {r1}")
    
    # check to see if the reference dir contains a .gtf and .fa file
    gtf_files = [f for f in os.listdir(reference_dir) if f.endswith('.gtf')]
    fa_files = [f for f in os.listdir(reference_dir) if f.endswith('.fa')]
    if len(gtf_files) != 1 or len(fa_files) != 1:
        raise Exception("The reference dir must have exactly one .gtf file and exactly one .fa file")
    
    # check to see if the output dir exists
    if not os.path.isdir(output_dir):
        raise Exception("Output directory does not exist")
    
    print("Inputs look OK")

    # build reference
    _build_rsem_reference(reference_dir)

def _trim_fastqs(input_dir: str, output_dir: str) -> None:
    os.makedirs(output_dir)
    fastqs = [f for f in os.listdir(input_dir) if '_R1_' in f and f.endswith('.fastq')]
    fastqs.sort()

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
        json_path = os.path.join(output_dir, html.replace('.html', '.json'))

        # run fastp
        os.system(f'fastp --in1 {fastq1_path} --in2 {fastq2_path} --out1 {fastq1_output_path} --out2 {fastq2_output_path} -h {html_path} -j {json_path}')
    
    # zip up the .html files for download
    # os.system(f'zip -r {output_dir}/__fastp_html_reports.zip {output_dir}/*.html')

def _summarize_fastp_reports(fastp_dir: str, output_dir: str, rsem_dir: str = "") -> None:
    if rsem_dir:
        aligned_counts = _read_rsem_counts(rsem_dir)

    os.makedirs(output_dir, exist_ok=True)
    fastp_reports = [f for f in os.listdir(fastp_dir) if f.endswith('.json')]
    fastp_reports.sort()

    read_counts_per_sample = dict()

    for report in fastp_reports:
        report_path = os.path.join(fastp_dir, report)

        # open the report and get the before/after filtering read counts
        with open(report_path, 'r') as f:
            report_data = json.load(f)
            total_reads_before = float(report_data['summary']['before_filtering']['total_reads'])
            total_reads_after = float(report_data['summary']['after_filtering']['total_reads'])
            low_quality_reads = float(report_data['filtering_result']['low_quality_reads'])
            
            # sample name is just the file name without the .json extension
            sample_name = report.replace('.json', '')

            read_counts_per_sample[sample_name] = {
                'total_reads_before_filtering': total_reads_before,
                'total_reads_after_filtering': total_reads_after,
                'low_quality_reads': low_quality_reads
            }

            # add the counts to the dictionary
            if rsem_dir:
                sample_name = report.replace('.json', '')
                
                if sample_name in aligned_counts:
                    read_counts_per_sample[sample_name]['aligned_reads'] = aligned_counts[sample_name]
                else:
                    print(f'Error: no RSEM counts found for {sample_name}')

    # write the summary to a .txt file
    if rsem_dir:
        summary_path = os.path.join(output_dir, '__fastp_rsem_summary.txt')
    else:
        summary_path = os.path.join(output_dir, '__fastp_summary.txt')

    progress_char = '█'
    with open(summary_path, 'w') as f:
        for sample, counts in read_counts_per_sample.items():
            f.write("Sample: " + sample + "\n")
            f.write("Total reads before filtering: " + str(counts['total_reads_before_filtering']) + "\n")
            f.write("Total reads after filtering: " + str(counts['total_reads_after_filtering']) + "\n")
            f.write("Low quality reads: " + str(counts['low_quality_reads']) + "\n")

            percentage_reads_kept = (counts['total_reads_after_filtering'] / counts['total_reads_before_filtering']) * 100
            percentage_reads_kept = round(percentage_reads_kept, 2)
            f.write("Percent of reads kept: " + str(percentage_reads_kept) + "%\n")

            if 'aligned_reads' in counts:
                f.write("Aligned reads: " + str(counts['aligned_reads']) + "\n")

                # calculate percentage rounded to 2 decimal places
                percentage_aligned = (counts['aligned_reads'] / counts['total_reads_after_filtering']) * 100
                percentage_aligned = round(percentage_aligned, 2)
                f.write("Percentage of reads aligned: " + str(percentage_aligned) + "%\n")

            # add a bar chart of the before/after filtering read counts
            before_filtering_bar_count = int(counts['total_reads_before_filtering'] / 2000000)
            after_filtering_bar_count = int(counts['total_reads_after_filtering'] / 2000000)

            before_bar = progress_char * before_filtering_bar_count
            after_bar = progress_char * after_filtering_bar_count

            # put a space every 10 ▉'s
            # i.e.,
            # ██████████ ██████████ ███
            before_bar = ' '.join([before_bar[i:i+10] for i in range(0, len(before_bar), 10)])
            after_bar = ' '.join([after_bar[i:i+10] for i in range(0, len(after_bar), 10)])

            f.write("Bar chart of reads before filtering: ---> " + before_bar + "\n")
            f.write("Bar chart of reads after filtering:  ---> " + after_bar + "\n")

            if 'aligned_reads' in counts:
                aligned_bar_count = int(counts['aligned_reads'] / 2000000)
                aligned_bar = progress_char * aligned_bar_count
                aligned_bar = ' '.join([aligned_bar[i:i+10] for i in range(0, len(aligned_bar), 10)])
                f.write("Bar chart of reads aligned:          ---> " + aligned_bar + "\n")

            f.write("\n\n")

def _run_rsem(input_dir: str, output_dir: str, genome_dir: str) -> None:
    # build star/rsem reference if it doesn't exist
    rsem_reference = _build_rsem_reference(genome_dir)

    os.makedirs(output_dir)

    # change working directory... fixes super annoying rsem bug.
    # RSEM tries to create a folder called {sample_name}.stat in the current working directory.
    # but in a docker container the cwd is not writable by default.
    # so we change the cwd to the output directory, which is presumably writable.
    # https://github.com/deweylab/RSEM/blob/8bc1e2115493c0cdf3c6bee80ef7a21a91b2acce/rsem-calculate-expression#L372
    original_wd = os.getcwd()
    os.chdir(output_dir)

    fastqs = [f for f in os.listdir(input_dir) if '_R1_' in f and f.endswith('.fastq')]
    fastqs.sort()

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

def _build_rsem_reference(genome_dir: str) -> str:
    gtf_files = [os.path.join(genome_dir, file) for file in os.listdir(genome_dir) if file.endswith(".gtf")]
    fa_files = [os.path.join(genome_dir, file) for file in os.listdir(genome_dir) if file.endswith(".fa")]

    star_reference_dir = os.path.join(genome_dir, "star_reference")
    star_reference_files = os.path.join(star_reference_dir, "star_reference")

    if not os.path.isdir(star_reference_dir):
        print("Building STAR reference")
        os.makedirs(star_reference_dir)
        os.system(f'rsem-prepare-reference --star -p 20 --gtf {gtf_files[0]} {fa_files[0]} {star_reference_files}')
    else:
        # TODO: check if the reference was built successfully by reading the Log.out file
        # if not, delete the directory and rebuild the reference

        # # read 'Log.out' file to see if the reference was built successfully
        # log_file = os.path.join(star_reference_dir, "Log.out")
        # with open(log_file, 'r') as f:
        #     lines = f.readlines()
        
        # # check for 'DONE: Genome generation, EXITING' in the log file
        # if not any('DONE: Genome generation, EXITING' in line for line in lines):
        #     # if the log file doesn't contain the 'DONE' message, then the reference wasn't built successfully
        #     raise Exception("Error building STAR reference")

        print("Found STAR reference directory. Skipping building reference.")
    
    return star_reference_files

def _write_counts_csv(rsem_dir: str, output_dir: str) -> str:
    os.makedirs(output_dir)
    # build the counts.csv file from the output in the RSEM directory
    os.system(f'Rscript /src/tximport.R --rsem_dir {rsem_dir} --output_dir {output_dir}')

def _read_rsem_counts(rsem_dir: str) -> dict:
    rsem_alignment_counts = dict()

    # get all dirs that end in .stat
    stat_dirs = [f for f in os.listdir(rsem_dir) if f.endswith('.stat')]
    stat_dirs.sort()
    for stat_dir in stat_dirs:
        stat_dir_path = os.path.join(rsem_dir, stat_dir)

        # sample name is just removing the .stat extension
        sample_name = stat_dir.replace('.stat', '')

        # read the .cnt file in the .stat directory
        cnt_file = [f for f in os.listdir(stat_dir_path) if f.endswith('.cnt')]

        if not cnt_file:
            print(f'Error: no .cnt file found in {stat_dir_path}')
            continue

        cnt_file_path = os.path.join(stat_dir_path, cnt_file[0])

        with open(cnt_file_path, 'r') as f:
            line0 = f.readline()
            # split by space
            split_line0 = line0.split()
            unalignable = int(split_line0[0]) * 2   # multiply by 2 because paired end
            alignable = int(split_line0[1]) * 2     # multiply by 2 because paired end
            rsem_alignment_counts[sample_name] = alignable
    
    return rsem_alignment_counts

def _make_output_dir(output_dir: str) -> str:
    current_datetime = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    output_dir = os.path.join(output_dir, current_datetime + "___rnatools")
    os.makedirs(output_dir)
    return output_dir

if __name__ == '__main__':
    main()