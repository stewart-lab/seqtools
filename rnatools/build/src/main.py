import os
import argparse
import json
from datetime import datetime

# TODO: make contrasts and such for DESeq2?

parser = argparse.ArgumentParser(description='GUPPy - Paired-end bulk RNA-Seq pipeline')
parser.add_argument('-i', '--fastq_dir', type=str, help='Directory containing input .fastq files.', required=True)
parser.add_argument('-r', '--reference_genome_dir', type=str, help='Directory containing the .gtf and .fa genome reference files.', required=True)
parser.add_argument('-o', '--output_dir', type=str, help='Output directory. By default this is the .fastq directory. This directory must exist.', required=False)
parser.add_argument('--resume', action='store_true', help='Resume a previous run of this pipeline. All args must be identical to resume successfully.', required=False)
args = parser.parse_args()

checkpoint_file = 'checkpoint.txt'
params_file = 'params.txt'

##############################################################################
## Runs the pipeline
##############################################################################
def main():
    fastq_dir = args.fastq_dir
    output_dir = args.output_dir
    genome_dir = args.reference_genome_dir
    resume = args.resume

    # set the output directory to the fastq directory if not specified
    if not output_dir:
        output_dir = fastq_dir

    # check for valid inputs before we start any processing
    _check_inputs(fastq_dir, genome_dir, output_dir)
    print("Inputs are valid!")

    # check if we are resuming a previous (interrupted) run
    if resume:
        print("Looking for previous run to resume...")
        existing_dir = _check_resume(output_dir, fastq_dir, genome_dir)
    else:
        existing_dir = None

    if not existing_dir:
        # create the timestamped output directory
        timestamped_outdir = _create_timestamped_outdir(output_dir)
        _write_checkpoint("STARTING PIPELINE", timestamped_outdir)
        _write_params(fastq_dir, genome_dir, timestamped_outdir)
    else:
        # resume from the existing directory
        print(f"Resuming run from {existing_dir}")
        timestamped_outdir = existing_dir
        _write_checkpoint("RESUMING PIPELINE", timestamped_outdir)

    # build STAR/RSEM reference if needed
    _run_rsem_prepare_reference(genome_dir, timestamped_outdir)

    # trim input fastqs
    fastp_dir = os.path.join(timestamped_outdir, '001_fastq_trimmed')
    _run_fastp(fastq_dir, fastp_dir, timestamped_outdir)

    # run rsem on trimmed fastqs
    rsem_dir = os.path.join(timestamped_outdir, '002_rsem')
    _run_rsem_calculate_expression(fastp_dir, rsem_dir, genome_dir, timestamped_outdir)

    # get counts matrix from RSEM output
    counts_dir = os.path.join(timestamped_outdir, '003_count_matrix')
    _run_tximport(rsem_dir, counts_dir, timestamped_outdir)

    # write all complete checkpoint
    _write_checkpoint("ALL STEPS COMPLETE", timestamped_outdir)


##############################################################################
## Checks pipeline inputs for validity
##############################################################################
def _check_inputs(fastq_dir: str, genome_dir: str, output_dir: str) -> None:
    # check to see if the output dir exists
    if not os.path.isdir(output_dir):
        raise Exception("Output directory does not exist")

    if not os.path.isdir(fastq_dir):
        raise Exception("Fastq directory does not exist")
    
    if not os.path.isdir(genome_dir):
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
    
    # check to see if the reference dir contains one .gtf and one .fa file
    gtf_files = [f for f in os.listdir(genome_dir) if f.endswith('.gtf')]
    fa_files = [f for f in os.listdir(genome_dir) if f.endswith('.fa')]
    if len(gtf_files) != 1 or len(fa_files) != 1:
        raise Exception("The reference dir must have exactly one .gtf file and exactly one .fa file")

##############################################################################
## Trim and filter .fastq files using fastp
##############################################################################
def _run_fastp(fastq_dir: str, fastp_dir: str, timestamped_outdir: str) -> None:
    os.makedirs(fastp_dir, exist_ok=True)

    fastqs = [f for f in os.listdir(fastq_dir) if '_R1_' in f and f.endswith('.fastq')]
    fastqs.sort()

    # get the list of fastq files that have already been trimmed (in case we are resuming a run)
    with open(os.path.join(timestamped_outdir, checkpoint_file), 'r') as f:
        checkpoint_lines = f.readlines()
    
    already_trimmed = [line.split(":")[1].strip() for line in checkpoint_lines if 'TRIMMED FASTQ FILES:' in line]
    
    # if all fastq files have already been trimmed, then we can skip this step
    if set(already_trimmed) == set(fastqs):
        print("All fastq files have already been trimmed. Skipping.")
        return

    # run fastp on each pair of fastq files
    for fastq1 in fastqs:
        if fastq1 in already_trimmed:
            print(f"Skipping {fastq1} because it has already been trimmed")
            continue

        # get input paths
        fastq2 = fastq1.replace('_R1_', '_R2_')
        html = fastq1.replace('.fastq', '.html').replace('_R1_', '').replace('001.html', '.html')
        fastq1_path = os.path.join(fastq_dir, fastq1)
        fastq2_path = os.path.join(fastq_dir, fastq2)

        # get output paths
        fastq1_output_path = os.path.join(fastp_dir, fastq1)
        fastq2_output_path = os.path.join(fastp_dir, fastq2)
        html_path = os.path.join(fastp_dir, html)
        json_path = os.path.join(fastp_dir, html.replace('.html', '.json'))

        # run fastp
        os.system(f'fastp --in1 {fastq1_path} --in2 {fastq2_path} --out1 {fastq1_output_path} --out2 {fastq2_output_path} -h {html_path} -j {json_path}')

        # summarize fastp reports
        _summarize_fastp_reports(fastp_dir, timestamped_outdir)
        _write_checkpoint(f"TRIMMED FASTQ FILES: {fastq1}", timestamped_outdir)
        _write_checkpoint(f"TRIMMED FASTQ FILES: {fastq2}", timestamped_outdir)
    
    _write_checkpoint("FASTP TRIMMING COMPLETE", timestamped_outdir)

    # zip up the .html files for download
    # os.system(f'zip -r {output_dir}/__fastp_html_reports.zip {output_dir}/*.html')

##############################################################################
## Builds the STAR reference from .gtf + .fa files
##############################################################################
def _run_rsem_prepare_reference(genome_dir: str, timestamped_outdir: str) -> str:
    gtf_files = [os.path.join(genome_dir, file) for file in os.listdir(genome_dir) if file.endswith(".gtf")]
    fa_files = [os.path.join(genome_dir, file) for file in os.listdir(genome_dir) if file.endswith(".fa")]

    star_reference_dir = os.path.join(genome_dir, "star_reference")
    star_reference_files = os.path.join(star_reference_dir, "star_reference")

    if not os.path.isdir(star_reference_dir):
        print("Building STAR reference")
        _write_checkpoint("BUILDING REFERENCE", timestamped_outdir)
        os.makedirs(star_reference_dir)
        os.system(f'rsem-prepare-reference --star -p 20 --gtf {gtf_files[0]} {fa_files[0]} {star_reference_files}')
        _write_checkpoint("REFERENCE BUILT", timestamped_outdir)
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
        _write_checkpoint("REFERENCE ALREADY EXISTS", timestamped_outdir)
    
    _write_checkpoint("REFERENCE OK", timestamped_outdir)
    return star_reference_files

##############################################################################
## Aligns the trimmed .fastq files to the STAR reference
##############################################################################
def _run_rsem_calculate_expression(fastq_dir: str, rsem_dir: str, genome_dir: str, timestamped_outdir: str) -> None:
    # get the STAR reference
    rsem_reference = _run_rsem_prepare_reference(genome_dir, timestamped_outdir)

    os.makedirs(rsem_dir, exist_ok=True)

    # change working directory... workaround for an RSEM bug.
    # RSEM tries to create a folder called {sample_name}.stat in the current working directory.
    # but in a docker container the cwd is not writable by default.
    # so we change the cwd to the output directory, which is presumably writable.
    # https://github.com/deweylab/RSEM/blob/8bc1e2115493c0cdf3c6bee80ef7a21a91b2acce/rsem-calculate-expression#L372
    original_wd = os.getcwd()
    os.chdir(rsem_dir)

    fastqs = [f for f in os.listdir(fastq_dir) if '_R1_' in f and f.endswith('.fastq')]
    fastqs.sort()

    # get the list of fastq files that have already been aligned (in case we are resuming a run)
    with open(os.path.join(timestamped_outdir, checkpoint_file), 'r') as f:
        checkpoint_lines = f.readlines()
    
    already_aligned = [line.split(":")[1].strip() for line in checkpoint_lines if 'RSEM ALIGNED:' in line]
    
    # if all fastq files have already been aligned, then we can skip this step
    if set(already_aligned) == set(fastqs):
        print("All fastq files have already been aligned. Skipping.")
        return

    rsem_temporary_dir = os.path.join(rsem_dir, 'rsem_temp')
    for fastq1 in fastqs:
        fastq2 = fastq1.replace('_R1_', '_R2_')
        sample_name = fastq1.split('_R1_')[0]

        if fastq1 in already_aligned:
            print(f"Skipping {fastq1} because it has already been aligned")
            continue

        fastq1_path = os.path.join(fastq_dir, fastq1)
        fastq2_path = os.path.join(fastq_dir, fastq2)
        os.system(f'rsem-calculate-expression --num-threads 40 --star --temporary-folder {rsem_temporary_dir} --paired-end {fastq1_path} {fastq2_path} {rsem_reference} {sample_name}')
        _write_checkpoint(f"RSEM ALIGNED: {fastq1}", timestamped_outdir)
        _write_checkpoint(f"RSEM ALIGNED: {fastq2}", timestamped_outdir)

        # summarize fastp reports again, this time with RSEM counts
        _summarize_fastp_reports(fastq_dir, timestamped_outdir, rsem_dir)

    # change working directory back to original working directory.
    # probably not necessary, but still.
    os.chdir(original_wd)

    _write_checkpoint("RSEM ALIGNMENT COMPLETE", timestamped_outdir)

##############################################################################
## Creates a counts.csv file from the RSEM output using tximport (R script)
##############################################################################
def _run_tximport(rsem_dir: str, counts_dir: str, timestamped_outdir: str) -> str:
    os.makedirs(counts_dir, exist_ok=True)
    # build the counts.csv file from the output in the RSEM directory
    os.system(f'Rscript /src/tximport.R --rsem_dir {rsem_dir} --output_dir {counts_dir}')
    _write_checkpoint("COUNTS MATRIX CREATED", timestamped_outdir)

##############################################################################
## Checks if there is a previous pipeline run to resume
##############################################################################
def _check_resume(output_dir: str, fastq_dir: str, genome_dir: str) -> str:
    # get all folders in the output directory
    all_files = os.listdir(output_dir)
    folders = [f for f in all_files if os.path.isdir(os.path.join(output_dir, f))]
    print(f"Found folders: {folders}")

    # reverse sort the folders so the most recent run is first
    folders.sort(reverse=True)

    # check if any of the folders contain the 'completed.txt' file
    for timestamped_outdir in folders:
        print(f"Checking {timestamped_outdir}")
        completed_file = os.path.join(output_dir, timestamped_outdir, checkpoint_file)
        params_file_path = os.path.join(output_dir, timestamped_outdir, params_file)
        same_params = False
        is_complete = False

        # check to see if the params has the same parameters as the current run
        current_run_params = _create_params_text(fastq_dir, genome_dir)
        if os.path.isfile(params_file_path):
            with open(params_file_path, 'r') as f:
                params_file_text = f.read()
            params_file_lines = [line.strip() for line in params_file_text.split('\n')]
            current_run_params_lines = [line.strip() for line in current_run_params.split('\n')]

            for i in range(len(current_run_params_lines)):
                if current_run_params_lines[i] == params_file_lines[i]:
                    same_params = True
                else:
                    print(f"Line {i + 1}")
                    print(f"Params lines don't match: {current_run_params_lines[i]} != {params_file_lines[i]}")
                    same_params = False
                    break

        # check to see if the completed_file has "ALL STEPS COMPLETE" in it
        if os.path.isfile(completed_file):
            with open(completed_file, 'r') as f:
                completed_text = f.read()
                if "ALL STEPS COMPLETE" in completed_text:
                    is_complete = True

        print(f"Same params: {same_params}")
        print(f"Is complete: {is_complete}")

        if same_params and not is_complete:
            return os.path.join(output_dir, timestamped_outdir)
        
    print("No previous run found to resume.")
    return None

##############################################################################
## Creates a summary of the fastp reports
##############################################################################
def _summarize_fastp_reports(fastp_dir: str, timestamped_outdir: str, rsem_dir: str = "") -> None:
    if rsem_dir:
        aligned_counts = _read_rsem_counts(rsem_dir)

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
                    print(f'WARNING: no RSEM counts found for {sample_name}')

    # write the summary to a .txt file
    if rsem_dir:
        summary_path = os.path.join(timestamped_outdir, '__fastp_rsem_summary.txt')
    else:
        summary_path = os.path.join(timestamped_outdir, '__fastp_summary.txt')

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

            # put a space every 10 █'s
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

##############################################################################
## Reads alignment rates from RSEM to add to the fastp report summary
##############################################################################
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
            print(f'WARNING: no .cnt file found in {stat_dir_path}')
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

##############################################################################
## Writes a params.txt file with information about this pipeline run
##############################################################################
def _write_params(fastq_dir: str, genome_dir: str, timestamped_outdir: str) -> None:
    # write a params file to the output directory
    params_path = os.path.join(timestamped_outdir, params_file)
    params_text = _create_params_text(fastq_dir, genome_dir)

    with open(params_path, 'w') as f:
        f.write(params_text)

        f.write("\n\nOutput directory:\n")
        f.write(timestamped_outdir)
    
    _write_checkpoint("WROTE PARAMS FILE", timestamped_outdir)

##############################################################################
## Summarizes information about this pipeline run for the params.txt file
##############################################################################
def _create_params_text(fastq_dir: str, genome_dir: str) -> None:
    # create text for a params file to write to the output directory
    # contains:
    # list of input .fastq files,
    # list of reference genome files,
    # list of tools used and their versions,
    # timestamped output directory,
    # this script
    conda_packages = os.popen('conda list').read()
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith('.gtf') or f.endswith('.fa')]
    genome_files.sort()
    fastq_files = [f for f in os.listdir(fastq_dir) if f.endswith('.fastq')]
    fastq_files.sort()

    pipeline_script_path = __file__
    with open(pipeline_script_path, 'r') as f:
        pipeline_script = f.read()

    # read the first few lines of the genome .gtf file to get the genome version
    fa_file = [f for f in genome_files if f.endswith('.gtf')][0]
    genome_version_lines = []
    with open(os.path.join(genome_dir, fa_file), 'r') as f:
        # read any line that begins with '#', up to 10 lines
        for i in range(10):
            line = f.readline()
            if line.startswith('#'):
                genome_version = line.strip()
                genome_version_lines.append(genome_version)


    # put the args in string form w/ json
    str_args = json.dumps(args.__dict__, indent=4)

    params_txt = f"""
    Conda packages:
    {conda_packages}
    \n\nReference genome files:
    {genome_files}
    \n\nGenome version:
    {genome_version_lines}
    \n\nFastq files:
    {fastq_files}
    \n\nArgs:
    {str_args}
    \n\nPipeline script:
    {pipeline_script}
    """

    return params_txt

##############################################################################
## Writes a message to the checkpoint.txt file
##############################################################################
def _write_checkpoint(message: str, timestamped_outdir: str) -> None:
    with open(os.path.join(timestamped_outdir, checkpoint_file), 'a') as f:
        f.write(message + '\n')

##############################################################################
## Makes a timestamped output directory
##############################################################################
def _create_timestamped_outdir(output_dir: str) -> str:
    current_datetime = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    timestamped_outdir = os.path.join(output_dir, current_datetime + "___rnatools")
    os.makedirs(timestamped_outdir)
    return timestamped_outdir

if __name__ == '__main__':
    main()