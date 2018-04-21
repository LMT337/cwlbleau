import os, csv, subprocess, datetime, argparse

# argument input woid
desc_str = """
        Program to parse cwl metrics.
    """
parser = argparse.ArgumentParser(description=desc_str)
parser.add_argument("-w", type=str, help='woid')
args = parser.parse_args()

# assign woid
if not args.w:
    print('usage:\n~awagner/bin/python3 cwlbleau.py -w <woid>')
    quit()
else:
    woid = args.w
    model_groups_id = 'model_groups.project.id=' + woid

# set working dir, results dic, date
working_dir = os.getcwd()
results = {}
results['WOID'] = woid
mm_dd_yy = datetime.datetime.now().strftime("%m%d%y")

# Fucntions pull metrics from last succeeded build dir
def verify_bamid(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.DictReader(infilecsv, delimiter='\t')
        for line in infile_reader:
            results['SEQ_ID'] = line['#SEQ_ID']
            results['FREEMIX'] = line['FREEMIX']

    return results


def insert_size_metrics(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        for line in infile_reader:
            if 'MEDIAN_INSERT_SIZE' in line:
                data = (next(infile_reader))
                results['MEDIAN_INSERT_SIZE'] = data[4]
                results['STANDARD_DEVIATION'] = data[5]

    return results


def flagstat_out(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv)
        mapped_rate = float()
        properly_paired_rate = float()
        for line in infile_reader:

            if 'mapped' in line[0] and '%' in line[0]:
                mapped_line = line[0]
                mapped_split = mapped_line.split('(')
                mapped_rate = float(mapped_split[1].split(':')[0].strip()[:-1])
                results['mapped_rate'] = mapped_rate

            if 'properly paired' in line[0]:
                properly_paired_line = line[0]
                properly_paired_split = properly_paired_line.split('(')
                properly_paired_rate = float(properly_paired_split[1].split(':')[0].strip()[:-1])
                results['properly_paired-rate'] = properly_paired_rate

            discordant_rate = mapped_rate - properly_paired_rate
            results['discordant_rate'] = discordant_rate

    return results


def mark_dups_metrics(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        percent_duplication = float()
        for line in infile_reader:
            if 'PERCENT_DUPLICATION' in line:
                data = next(infile_reader)
                percent_duplication = float(data[8])
                results['PERCENT_DUPLICATION'] = data[8]

    return percent_duplication


def gcbias_metrics_summary(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        for line in infile_reader:
            if 'ALIGNED_READS' in line:
                data = next(infile_reader)
                results['ALIGNED_READS'] = data[4]

    return results


def alignment_summary_metrics(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        pf_aligned_bases = float()
        for line in infile_reader:
            if 'FIRST_OF_PAIR' in line:
                results['FOP: PF_MISMATCH_RATE'] = line[12]
            if 'SECOND_OF_PAIR' in line:
                results['SOP: PF_MISMATCH_RATE'] = line[12]
            if 'PAIR' in line and not '_' in line:
                pf_aligned_bases = int(line[7])
                results['TOTAL_READS'] = line[1]
                results['PF_READS'] = line[2]
                results['PF_READS_ALIGNED'] = line[5]
                results['PF_ALIGNED_BASES'] = line[7]
                results['PF_HQ_ALIGNED_Q20_BASE'] = line[10]
                results['PCT_ADAPTER'] = line[23]

    return pf_aligned_bases


def wgs_metrics(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        genome_territory = int()
        for line in infile_reader:
            if 'GENOME_TERRITORY' in line:
                data = next(infile_reader)
                genome_territory = int(data[0])
                results['GENOME_TERRITORY'] = data[0]
                results['MEAN_COVERAGE'] = data[1]
                results['SD_COVERAGE'] = data[2]
                results['PCT_10X'] = data[14]
                results['PCT_20X'] = data[16]
                results['PCT_30X'] = data[18]
                results['HET_SNP_SENSITIVITY'] = data[26]
                results['HET_SNP_Q'] = data[27]

    return genome_territory


# run genome model command to generate model info
model_groups = subprocess.check_output(['genome', 'model', 'list', model_groups_id, "--show",
                                        "last_succeeded_build.id,name,status,last_succeeded_build.data_directory",
                                        "--style=tsv", "--nohead"]).decode('utf-8').splitlines()

# header for outfile
metrics_header = ['WOID', 'last_succeeded_build', 'name', 'status', 'data_directory', 'properly_paired-rate','PF_READS',
                  'FREEMIX', 'discordant_rate', 'FOP: PF_MISMATCH_RATE', 'GENOME_TERRITORY', 'mapped_rate',
                  'SD_COVERAGE', 'HAPLOID COVERAGE', 'TOTAL_READS', 'PF_READS_ALIGNED', 'SEQ_ID', 'HET_SNP_SENSITIVITY',
                  'MEDIAN_INSERT_SIZE', 'PCT_20X', 'PF_ALIGNED_BASES', 'PCT_30X', 'PERCENT_DUPLICATION', 'PCT_ADAPTER',
                  'ALIGNED_READS', 'PCT_10X', 'STANDARD_DEVIATION', 'MEAN_COVERAGE', 'PF_HQ_ALIGNED_Q20_BASE',
                  'SOP: PF_MISMATCH_RATE', 'HET_SNP_Q']

# outfile
cwd_metrics_outfile = woid + '.cwl.metrics.' + mm_dd_yy + '.tsv'

# call methods to generate results
with open(cwd_metrics_outfile, 'w') as outfilecsv:

    metrics_writer = csv.DictWriter(outfilecsv, fieldnames=metrics_header, delimiter='\t')
    metrics_writer.writeheader()

    for line in model_groups:

        info = line.split('\t')

        if 'Succeeded' in info[2]:
            results['last_succeeded_build'] = info[0]
            results['name'] = info[1]
            results['status'] = info[2]
            results['data_directory'] = info[3]

            os.chdir(info[3] + '/results')

            verify_bamid('VerifyBamId.selfSM')
            insert_size_metrics('InsertSizeMetrics.txt')
            flagstat_out('flagstat.out')
            perc_dup = mark_dups_metrics('mark_dups_metrics.txt')
            gcbias_metrics_summary('GcBiasMetricsSummary.txt')
            pfalgnbases = alignment_summary_metrics('AlignmentSummaryMetrics.txt')
            genome_terr = wgs_metrics('WgsMetrics.txt')

            haploid_coverage = pfalgnbases * ((1 - perc_dup)/genome_terr)
            results['HAPLOID COVERAGE'] = haploid_coverage

            os.chdir(working_dir)

            metrics_writer.writerow(results)

        else:

            print('{} {} Build Failed'.format(info[0], info[1]))