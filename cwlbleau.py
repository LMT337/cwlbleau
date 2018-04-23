import os, csv, subprocess, datetime, argparse

# argument input woid
desc_str = """
        Program to parse cwl metrics.
    """
parser = argparse.ArgumentParser(description=desc_str)
parser.add_argument("-w", type=str, help='woid')
parser.add_argument("-f", type=str, help='file of woid\'s (without header)')
args = parser.parse_args()

woid_list = []
# assign woid
if not args.w and not args.f:
    print('usage:\n~awagner/bin/python3 cwlbleau.py -w <woid>\n~awagner/bin/python3 cwlbleau.py -f <file of woids>')
    quit()
elif args.w:
    woid = args.w
    woid_list.append(woid)
    model_groups_id = 'model_groups.project.id=' + woid
elif args.f:
    with open(args.f, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv)
        for line in infile_reader:
            woid_list.append(line[0])

# set working dir, results dic, date
working_dir = os.getcwd()
results = {}
mm_dd_yy = datetime.datetime.now().strftime("%m%d%y")
mmddyy_slash = datetime.datetime.now().strftime("%m/%d/%y")


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
            if 'MEAN_INSERT_SIZE' in line:
                data = (next(infile_reader))
                results['MEAN_INSERT_SIZE'] = data[4]
                results['STANDARD_DEVIATION'] = data[5]

    return results


def flagstat_out(infile):
    with open(infile, 'r') as infilecsv:

        infile_reader = csv.reader(infilecsv)

        mapped_rate = float()
        mapped_int = float()
        properly_paired_rate = float()
        wmmtadc = float()

        for line in infile_reader:

            if 'mapped' in line[0] and '% : N/A' in line[0]:
                mapped_line = line[0]
                mapped_split = mapped_line.split('(')
                mapped_rate = float(mapped_split[1].split(':')[0].strip()[:-1])
                mapped_int = float(mapped_split[0].split('+')[0].strip())
                results['mapped_rate'] = mapped_rate

            if 'properly paired' in line[0]:
                properly_paired_line = line[0]
                properly_paired_split = properly_paired_line.split('(')
                properly_paired_rate = float(properly_paired_split[1].split(':')[0].strip()[:-1])
                results['properly_paired-rate'] = properly_paired_rate

            if 'with mate mapped to a different chr' in line[0] and 'mapQ>=' not in line[0]:
                wmmtadc_line = line[0]
                wmmtadc = float(wmmtadc_line.split('+')[0].strip())

        discordant_rate = mapped_rate - properly_paired_rate
        results['discordant_rate'] = discordant_rate

        inter_chromosomal_pairing_rate = wmmtadc / mapped_int
        results['inter-chromosomal_Pairing rate'] = inter_chromosomal_pairing_rate

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
                                        "last_succeeded_build.id,name,status,last_succeeded_build.data_directory,"
                                        "subject.name", "--style=tsv", "--nohead"]).decode('utf-8').splitlines()


metrics_header = ['Admin', 'WorkOrder','date_QC','sample_name','model_name','last_succeeded_build','data_directory',
                  'cram_file','status', 'ALIGNED_READS','mapped_rate','FOP: PF_MISMATCH_RATE','SOP: PF_MISMATCH_RATE',
                  'FREEMIX','HAPLOID COVERAGE','PCT_10X', 'PCT_20X','PCT_30X','discordant_rate',
                  'inter-chromosomal_Pairing rate','HET_SNP_Q','HET_SNP_SENSITIVITY', 'HET_SNP_SENSITIVITY',
                  'MEAN_COVERAGE','SD_COVERAGE','MEAN_INSERT_SIZE','STANDARD_DEVIATION','PCT_ADAPTER','PF_READS',
                  'PF_ALIGNED_BASES','PERCENT_DUPLICATION','TOTAL_READS','properly_paired-rate',
                  'PF_HQ_ALIGNED_Q20_BASE','PF_READS_ALIGNED','GENOME_TERRITORY','SEQ_ID']

for woid in woid_list:
    cwd_metrics_outfile = woid + '.cwl.metrics.' + mm_dd_yy + '.tsv'

    #Admin name
    admin_collections = subprocess.check_output(["wo_info", "--report", "billing", "--woid", woid]).decode(
        'utf-8').splitlines()
    for ap in admin_collections:
        if 'Administration Project' in ap:
            ap_new = ap.split(':')[1].strip()
            results['Admin'] = ap_new

    # call methods to generate results
    with open(cwd_metrics_outfile, 'w') as outfilecsv:

        metrics_writer = csv.DictWriter(outfilecsv, fieldnames=metrics_header, delimiter='\t')
        metrics_writer.writeheader()

        for line in model_groups:

            info = line.split('\t')

            if 'Succeeded' in info[2]:

                results['WorkOrder'] = woid
                results['date_QC'] = mmddyy_slash
                results['last_succeeded_build'] = info[0]
                results['model_name'] = info[1]
                results['status'] = info[2]
                results['data_directory'] = info[3]
                results['sample_name'] = info[4]

                os.chdir(info[3] + '/results')

                results['cram_file'] = 'NA'
                if os.path.isfile('final.cram'):
                    results['cram_file'] = os.getcwd() + '/final.cram'

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