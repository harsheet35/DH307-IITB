GENOMES_DIR='/home/user/Desktop/Prakhar/genome'
OUT_DIR = '/home/user/Desktop/Prakhar/my_rnaseq_analysis/rna_seq'
RAWDATA_DIR ='/home/user/Desktop/Prakhar/my_rnaseq_analysis/raw_data'
SAMPLES=['SRR1672666']#, 'SRR1672666'

GENOME_BUILD = 'GRCh38'
GENOME_FASTA = GENOMES_DIR + '/' + 'genome.fa'
SRC_DIR = '/home/user/Desktop/Prakhar/my_rnaseq_analysis/scripts'
STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/star_annotated'
GTF = GENOMES_DIR + '/' 'genome_ann.gtf'

GENE_NAMES = '/home/user/Desktop/Prakhar/my_rnaseq_analysis/gene_names.tsv'
GENE_LENGTHS = '/home/user/Desktop/Prakhar/my_rnaseq_analysis/gene_lengths.tsv'  #+ GENOME_BUILD+'gene_lengths.tsv'
GENE_NAME_MAP = '/home/user/Desktop/Prakhar/my_rnaseq_analysis/gene_name_map.tsv'
