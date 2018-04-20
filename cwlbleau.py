import os, csv

os.chdir('results')
results = {}


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


verify_bamid('VerifyBamId.selfSM')
insert_size_metrics('InsertSizeMetrics.txt')
flagstat_out('flagstat.out')
perc_dup = mark_dups_metrics('mark_dups_metrics.txt')
gcbias_metrics_summary('GcBiasMetricsSummary.txt')
pfalgnbases = alignment_summary_metrics('AlignmentSummaryMetrics.txt')
genome_terr = wgs_metrics('WgsMetrics.txt')

haploid_coverage = pfalgnbases * ((1 - perc_dup)/genome_terr)
results['HAPLOID COVERAGE'] = haploid_coverage