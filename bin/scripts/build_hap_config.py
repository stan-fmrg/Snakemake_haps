import argparse
import os
import subprocess
#import yaml




# Read in some command line arguments

parser = argparse.ArgumentParser()
parser.add_argument("-se", "--se_samples", type=str, default = "data/config/se_samples.csv", help="Single-end sample file, single_end.csv")
parser.add_argument("-pe", "--pe_samples", type=str, default = "data/config/pe_samples.csv", help="Paired-end sample file, paired_end.csv")
parser.add_argument("-n", "--run_name", type=str, default = "haps_run", help="Name of run")
parser.add_argument("-r", "--scripts_dir", type=str, default = "bin/scripts", help="Relative path for scripts")
parser.add_argument("-g", "--genome", type=str, default = "data/genome/sample.fasta", help="Relative path for reference genome")
parser.add_argument("-c", "--config_dir", type=str, default = "data/config", help="Relative path for various config files")
parser.add_argument("-d", "--data_dir", type=str, default = "data/fastq_raw", help="Relative path for raw fastq files")
parser.add_argument("-a", "--adapter_dir", type=str, default = "data/adapters", help="Relative path for adapter files")
parser.add_argument("-s", "--split_genome", type=int, default = 5, help="Split the reference genome into n pieces for parallel variant calling")
parser.add_argument("-v", "--diff_vcf", type=str, default = "", help="Assign a vcf file with known variants for comparison with the new vcf for common SNPs")
parser.add_argument("-p", "--run_haplotyper", type=str, default = True, help="Allow rad_haplotyper to design haplotypes from VCF files (Default value = True)")
parser.add_argument("-q", "--Quality_control", type=str, default= True, help="Allow FastQC tool to proceed with the quality control of the raw reads and the quality trimmed reads (Default value = True)")
parser.add_argument("-f", "--run_FLASH", type=str, default = True, help="FLASH merges paired-end reads into a single one. FLASH is designed to merge pairs of reads when the original DNA fragments are shorter than twice the length of reads. (Default value = True)")
args = parser.parse_args()


subprocess.call(['bash','./bin/scripts/Sample_list.sh',args.data_dir])

# Function for autovivification of nested python dictionaries:
# https://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries?noredirect=1&lq=1
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

# Create empty dicts for storing the data
samps = Vividict()
yaml_dict = Vividict()

yaml_dict['run_name'] = args.run_name
yaml_dict['scripts_dir'] = args.scripts_dir
yaml_dict['config_dir'] = args.config_dir
yaml_dict['genome'] = args.genome
yaml_dict['data_dir'] = args.data_dir
yaml_dict['adapter_dir'] = args.adapter_dir
yaml_dict['home_dir'] = os.getcwd()


# Parse the sample file and build data structures
if (args.pe_samples):
    with open(args.pe_samples) as SAMP:
        samp_header = SAMP.readline()
        for line in SAMP:
            fields = line.rstrip().split(",")
            yaml_dict['samples_pe'][fields[0]]['forward'] = args.data_dir + "/" + str(fields[1])
            yaml_dict['samples_pe'][fields[0]]['reverse'] = args.data_dir + "/" + str(fields[2])
            yaml_dict['samples_pe'][fields[0]]['adapter'] = args.adapter_dir + "/" + fields[3]
            yaml_dict['samples_pe'][fields[0]]['phred'] = fields[4]

if (args.se_samples):
    with open(args.se_samples) as SAMP:
        samp_header = SAMP.readline()
        for line in SAMP:
            fields = line.rstrip().split(",")
            yaml_dict['samples_se'][fields[0]]['forward'] = args.data_dir + "/" + str(fields[1])
            yaml_dict['samples_se'][fields[0]]['adapter'] = args.adapter_dir + "/" + fields[2]
            yaml_dict['samples_se'][fields[0]]['phred'] = fields[3]
        
#print samps

# Write the snakemake config file
with open("config.yaml", mode = "w") as OUT:
    OUT.write('{0!s} {1!s}\n'.format("run_name:", yaml_dict['run_name']))
    OUT.write('{0!s} {1!s}\n'.format("home_dir:", yaml_dict['home_dir']))
    OUT.write('{0!s} {1!s}\n'.format("scripts_dir:", yaml_dict['scripts_dir']))
    OUT.write('{0!s} {1!s}\n'.format("config_dir:", yaml_dict['config_dir']))
    OUT.write('{0!s} {1!s}\n'.format("genome:", yaml_dict['genome']))
    
    OUT.write('{0!s}\n'.format("###################################################\n# Trimmomatic read trimming\n# :Parameters:\n# - Trimmomatic: any parameters recognised by Trimmomatic tool"))
    OUT.write('{0!s}\n'.format("Read Trimming:"))
    OUT.write('  {0!s} {1!s}\n'.format("trimmomatic_params:", "LEADING:20 TRAILING:20 AVGQUAL:30 MINLEN:100 TOPHRED33 SLIDINGWINDOW:4:15\n"))
    
    OUT.write('{0!s}\n'.format("###################################################\n# Samtools Fixmate\n# :Parameters:\n# - Samtools Fixmate: any parameters recognised by Samtools Fixmate:"))
    OUT.write('{0!s}\n'.format("SamTools_Fixmate:"))
    OUT.write('  {0!s} {1!s}\n'.format("fixmate_params:", "-r -m\n"))
    
    OUT.write('{0!s}\n'.format("###################################################\n# Samtools MarkDup\n# :Parameters:\n# - Samtools MarkDup: any parameters recognised by Samtools MarkDup:"))
    OUT.write('{0!s}\n'.format("Marking_Duplicates:"))
    OUT.write('  {0!s} {1!s}\n'.format("markdup_params:", "-s\n"))
    
    OUT.write('{0!s}\n'.format("###################################################\n# BWA - Mapping\n# :Parameters:\n# - BWA_mem_params: any parameters recognised by BWA MEM tool"))
    OUT.write('{0!s}\n'.format("BWA-Mapping:"))
    OUT.write('  {0!s} {1!s}\n'.format("BWA_mem_params:", "-t 8, -B 4, -O 6\n"))
    
    OUT.write('{0!s}\n'.format("###################################################\n# Samtools Read Filtering - Mapping\n# :Parameters:\n# - Samtools view: any parameters recognised by Samtools view"))
    OUT.write('{0!s}\n'.format("Samtools_read_filtering:"))
    OUT.write('  {0!s} {1!s}\n'.format("samtools_view_filtering:", "-f 3 -F 780\n"))
    
    OUT.write('{0!s}\n'.format("###################################################\n# FreeBayes Variant Calling\n# :Parameters:\n# - FreeBayes: any parameters recognised by FreeBayes"))
    OUT.write('{0!s}\n'.format("FreeBayes_variant_calling:"))
    OUT.write('  {0!s} {1!s}\n'.format("freebayes_params:", "--haplotype-length 0 --no-complex --no-mnps --no-indels --use-duplicate-reads\n"))

    OUT.write('{0!s}\n'.format("###################################################\n# Variant Filtering \n# :Parameters:\n# - VCFtools: any parameters recognised by VCFtools"))
    OUT.write('{0!s}\n'.format("VCFtools:"))
    OUT.write('  {0!s} {1!s}\n'.format("variant_filtering:", "--minDP 10 --minGQ 20 --minQ 20 --max-missing 0.25 --maf 0.05 --recode --recode-INFO-all\n"))

    OUT.write('{0!s}\n'.format("###################################################\n# Samtools mapping parameters \n# :Parameters:\n# - Samtools view: any parameters recognised by Samtools view"))
    OUT.write('{0!s}\n'.format("Samtools_mapping_params:"))
    OUT.write('  {0!s} {1!s}\n'.format("samtools_view_params:", "-Sbhu -@ 8\n"))

    OUT.write('{0!s}\n'.format("###################################################\n# FLASH parameters \n# :Parameters:\n# - FLASH: any parameters recognised by FLASH"))
    OUT.write('{0!s}\n'.format("FLASH:"))
    OUT.write('  {0!s} {1!s}\n'.format("FLASH_params:", "-m 30 -M 250\n"))
    
    OUT.write('{0!s}\n'.format("###################################################\n# Haplotyping parameters \n# :Parameters:\n# - rad_haplotyper: any parameters recognised by rad_haplotyper"))
    OUT.write('{0!s}\n'.format("Haplotyping_params:"))
    OUT.write('  {0!s} {1!s}\n'.format("rad_haplotyper_params:", "-m 0.75 -z 3\n"))

    OUT.write('{0!s}\n'.format("###################################################\n# Using a pre-existing VCF file with variants to find overlapping SNPs from the snakemake VCF file and use these for haplotypes \n# :Parameters:\n# - Relative path of VCF file (Default: None)"))
    OUT.write('{0!s} {1!s}\n'.format("Validation:", args.diff_vcf))

    OUT.write('{0!s}\n'.format("###################################################\n# Quality control of fastq reads before and after quality trimming (Opional) \n# :Parameters:\n# - True or False (Default: True)"))
    OUT.write('{0!s} {1!s}\n'.format("QC_Tests:", args.Quality_control))

    OUT.write('{0!s}\n'.format("###################################################\n# Executing rad_haplotyper using the filtered VCF file from snakemake workflow (Optional) \n# :Parameters:\n# - True or False (Default: True)"))
    OUT.write('{0!s} {1!s}\n'.format("Haplotyper:", args.run_haplotyper))

    OUT.write('{0!s}\n'.format("###################################################\n# merges paired-end reads into a single one. FLASH is designed to merge pairs of reads when the original DNA fragments are shorter than twice the length of reads (Optional) \n# :Parameters:\n# - True or False (Default: True)"))
    OUT.write('{0!s} {1!s}\n'.format("Flash:", args.run_FLASH))

    OUT.write('{0!s}\n'.format("###################################################\n# Splitting reference into parts for multi-threading and faster analysis \n# :Parameters:\n# - No. of parts (Default: 5)"))
    OUT.write('{0!s}\n'.format("variant_call_params:"))
    OUT.write('  {0!s} {1!s}\n'.format("split_genome:", args.split_genome))
    OUT.write('{0!s}\n'.format("###################################################"))
    if (args.se_samples):
        OUT.write('{0!s}\n'.format("illumina_se:"))
        for sample in yaml_dict['samples_se']:
            OUT.write('  {0!s}{1!s}\n'.format(sample, ":"))
            OUT.write('    {0!s} {1!s}\n'.format("single:", yaml_dict['samples_se'][sample]['forward']))
            OUT.write('    {0!s} {1!s}\n'.format("phred:", yaml_dict['samples_se'][sample]['phred']))
            OUT.write('    {0!s} {1!s}\n'.format("adapter:", yaml_dict['samples_se'][sample]['adapter']))
    if (args.pe_samples):
        OUT.write('{0!s}\n'.format("illumina_pe:"))
        for sample in yaml_dict['samples_pe']:
            OUT.write('  {0!s}{1!s}\n'.format(sample, ":"))
            OUT.write('    {0!s} {1!s}\n'.format("forward:", yaml_dict['samples_pe'][sample]['forward']))
            OUT.write('    {0!s} {1!s}\n'.format("reverse:", yaml_dict['samples_pe'][sample]['reverse']))
            OUT.write('    {0!s} {1!s}\n'.format("phred:", yaml_dict['samples_pe'][sample]['phred']))
            OUT.write('    {0!s} {1!s}\n'.format("adapter:", yaml_dict['samples_pe'][sample]['adapter']))

OUT.close()
