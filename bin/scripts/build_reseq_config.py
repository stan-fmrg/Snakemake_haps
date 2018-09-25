import argparse
import os
#import yaml

# Read in some command line arguments

parser = argparse.ArgumentParser()
parser.add_argument("-se", "--se_samples", type=str, help="Single-end sample file")
parser.add_argument("-pe", "--pe_samples", type=str, help="Paired-end sample file")
parser.add_argument("-n", "--run_name", type=str, default = "reseq_run", help="Name of run")
parser.add_argument("-r", "--scripts_dir", type=str, default = "bin/scripts", help="Relative path for scripts")
parser.add_argument("-g", "--genome", type=str, default = "data/genome/sample.fasta", help="Relative path for reference genome")
parser.add_argument("-c", "--config_dir", type=str, default = "data/config", help="Relative path for various config files")
parser.add_argument("-d", "--data_dir", type=str, default = "data/fastq_raw", help="Relative path for raw fastq files")
parser.add_argument("-a", "--adapter_dir", type=str, default = "data/adapters", help="Relative path for adapter files")
args = parser.parse_args()

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
    OUT.write('{0!s} {1!s}\n'.format("trimmomatic_params:", "LEADING:20 TRAILING:20 AVGQUAL:30 MINLEN:32 TOPHRED33"))
    if (args.se_samples):
        OUT.write('{0!s}\n'.format("illumina_se:"))
        for sample in yaml_dict['samples_se']:
            OUT.write('  {0!s}{1!s}\n'.format(sample, ":"))
            OUT.write('    {0!s} {1!s}\n'.format("forward:", yaml_dict['samples_se'][sample]['forward']))
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
