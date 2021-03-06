import os, csv, subprocess, datetime, argparse

# argument input
desc_str = """
        Program to parse cwl metrics.
    """
parser = argparse.ArgumentParser(description=desc_str)
parser.add_argument("-w", type=str, help='woid')
parser.add_argument("-fw", type=str, help='file of woid\'s (without header)')
parser.add_argument("-a", type=str, help='anp')
parser.add_argument("-fa", type=str, help='file of AnP\'s (without header)')
args = parser.parse_args()

id_list = []

# assign woid
if not args.w and not args.a and not args.fw and not args.fa:
    print('usage:\n~awagner/bin/python3 cwlbleau.py -w <woid>\n~awagner/bin/python3 cwlbleau.py -f <file of woids>')
    quit()
elif args.w:
    anp_or_woid = "WorkOrder"
    id_list.append(args.w)
elif args.a:
    anp_or_woid = "AnP"
    id_list.append(args.a)
elif args.fw:
    anp_or_woid = "WorkOrder"
    with open(args.fw, 'r') as infilecsv:
        for line in infilecsv:
            id_list.append(line.rstrip())
elif args.fa:
    anp_or_woid = "AnP"
    with open(args.fa) as infilecsv:
        for line in infilecsv:
            id_list.append(line.rstrip())

# set working dir, results dic, date
working_dir = os.getcwd()
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
                results['MEAN_INSERT_SIZE'] = data[5]
                results['STANDARD_DEVIATION'] = data[6]
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


def hs_metrics(infile):
    with open(infile, 'r') as infilecsv:
        infile_reader = csv.reader(infilecsv, delimiter='\t')
        for line in infile_reader:
            if 'BAIT_SET' in line:
                hs_metrics_header = line
                # hs_metrics_data_one = next(infile_reader)
                hs_metrics_data_two = next(infile_reader)
                hs_metrics_dict = dict(zip(hs_metrics_header, hs_metrics_data_two))
        for metric in hs_metrics_dict:
            results[metric] = hs_metrics_dict[metric]
    return hs_metrics_header


def write_results(results_dict, outfile, header_list):
    if not os.path.isfile(outfile):
        with open(cwd_metrics_outfile, 'w') as outfilecsv:
            metrics_writer = csv.DictWriter(outfilecsv, fieldnames=header_list, delimiter='\t')
            metrics_writer.writeheader()
            metrics_writer.writerow(results_dict)
    elif os.path.isfile(outfile):
        with open(cwd_metrics_outfile, 'a') as outfilecsv:
            metrics_writer = csv.DictWriter(outfilecsv, fieldnames=header_list, delimiter='\t')
            metrics_writer.writerow(results_dict)
    return


met_wgs_header = ['Admin', 'WorkOrder','date_QC','sample_name','model_name','last_succeeded_build','data_directory',
                  'cram_file','status', 'ALIGNED_READS','mapped_rate','FOP: PF_MISMATCH_RATE','SOP: PF_MISMATCH_RATE',
                  'FREEMIX','HAPLOID COVERAGE','PCT_10X', 'PCT_20X','PCT_30X','discordant_rate',
                  'inter-chromosomal_Pairing rate','HET_SNP_Q','HET_SNP_SENSITIVITY',
                  'MEAN_COVERAGE','SD_COVERAGE','MEAN_INSERT_SIZE','STANDARD_DEVIATION','PCT_ADAPTER','PF_READS',
                  'PF_ALIGNED_BASES','PERCENT_DUPLICATION','TOTAL_READS','properly_paired-rate',
                  'PF_HQ_ALIGNED_Q20_BASE','PF_READS_ALIGNED','GENOME_TERRITORY','SEQ_ID']
met_wgs_header[1] = anp_or_woid


for id in id_list:

    if args.w or args.fw:
        print('cwlbleau\'ing: {}'.format(id))
        model_groups_id = 'model_groups.project.id=' + id

        # run genome model command to generate model info
        model_groups = subprocess.check_output(['genome', 'model', 'list', model_groups_id, "--show",
                                                "last_succeeded_build.id,name,status,last_succeeded_build.data_directory,"
                                                "subject.name", "--style=tsv", "--nohead"]).decode('utf-8').splitlines()

    if args.a or args.fa:
        print('cwlbleau\'ing: {}'.format(id))
        model_groups_id = 'analysis_project.id=' + id

        # run genome model command to generate model info
        model_groups = subprocess.check_output(['genome', 'model', 'list', model_groups_id, "--show",
                                                "last_succeeded_build.id,name,status,last_succeeded_build.data_directory,"
                                                "subject.name", "--style=tsv", "--nohead"]).decode('utf-8').splitlines()

    cwd_metrics_outfile = id + '.cwl.metrics.' + mm_dd_yy + '.tsv'
    print('Metrics outfile: {}'.format(cwd_metrics_outfile))
    if os.path.isfile(cwd_metrics_outfile):
        os.remove(cwd_metrics_outfile)

    # Admin name
    ap_new = "NA"
    if args.w or args.fw:
        admin_collections = subprocess.check_output(["wo_info", "--report", "billing", "--woid", id]).decode(
            'utf-8').splitlines()
        for ap in admin_collections:
            if 'Administration Project' in ap:
                ap_new = ap.split(':')[1].strip()

    # call methods to generate results
    for line in model_groups:

        results = {}

        info = line.split('\t')

        if 'Succeeded' in info[2]:

            results.clear()
            results['Admin'] = ap_new
            results[anp_or_woid] = id
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

            if os.path.isfile('VerifyBamId.selfSM'):
                verify_bamid('VerifyBamId.selfSM')
            else:
                results['SEQ_ID'] = 'FNF'
                results['FREEMIX'] = 'FNF'

            if os.path.isfile('InsertSizeMetrics.txt'):
                insert_size_metrics('InsertSizeMetrics.txt')
            else:
                results['MEAN_INSERT_SIZE'] = 'FNF'
                results['STANDARD_DEVIATION'] = 'FNF'

            if os.path.isfile('flagstat.out'):
                flagstat_out('flagstat.out')
            else:
                results['mapped_rate'] = 'FNF'
                results['properly_paired-rate'] = 'FNF'
                results['discordant_rate'] = 'FNF'
                results['inter-chromosomal_Pairing rate'] = 'FNF'

            if os.path.isfile('mark_dups_metrics.txt'):
                perc_dup = mark_dups_metrics('mark_dups_metrics.txt')
            else:
                results['PERCENT_DUPLICATION'] = 'FNF'
                perc_dup = False

            if os.path.isfile('GcBiasMetricsSummary.txt'):
                gcbias_metrics_summary('GcBiasMetricsSummary.txt')
            else:
                results['ALIGNED_READS'] = 'FNF'

            if os.path.isfile('AlignmentSummaryMetrics.txt'):
                pfalgnbases = alignment_summary_metrics('AlignmentSummaryMetrics.txt')
            else:
                results['TOTAL_READS'] = 'FNF'
                results['PF_READS'] = 'FNF'
                results['PF_READS_ALIGNED'] = 'FNF'
                results['PF_ALIGNED_BASES'] = 'FNF'
                results['PF_HQ_ALIGNED_Q20_BASE'] = 'FNF'
                results['PCT_ADAPTER'] = 'FNF'
                pfalgnbases = False

            if os.path.isfile('WgsMetrics.txt'):
                genome_terr = wgs_metrics('WgsMetrics.txt')
            else:
                results['GENOME_TERRITORY'] = 'FNF'
                results['MEAN_COVERAGE'] = 'FNF'
                results['SD_COVERAGE'] = 'FNF'
                results['PCT_10X'] = 'FNF'
                results['PCT_20X'] = 'FNF'
                results['PCT_30X'] = 'FNF'
                results['HET_SNP_SENSITIVITY'] = 'FNF'
                results['HET_SNP_Q'] = 'FNF'
                genome_terr = False

            if pfalgnbases and perc_dup and genome_terr:
                haploid_coverage = pfalgnbases * ((1 - perc_dup)/genome_terr)
                results['HAPLOID COVERAGE'] = haploid_coverage
            else:
                results['HAPLOID COVERAGE'] = 'FNF'

            if os.path.isfile('HsMetrics.txt'):
                header_add_fields = hs_metrics('HsMetrics.txt')
                metrics_header = met_wgs_header + header_add_fields
                # metrics_header = list(results.keys())
                # metrics_header.sort()
            else:
                metrics_header = met_wgs_header

            os.chdir(working_dir)
            write_results(results, cwd_metrics_outfile, metrics_header)

        else:
            print('{} {} Build Failed'.format(info[0], info[1]))