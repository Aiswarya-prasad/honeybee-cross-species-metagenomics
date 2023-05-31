import os
import sys
import argparse
import gzip
from Bio import SeqIO

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--flagstat', nargs="+", metavar="flagstat",required=True, help="space-separated list of paths to flagstat files", action="store")
requiredNamed.add_argument('--scaffolds', nargs="+", metavar="scaffolds",required=True, help="space-separated list of paths to scaffolds files", action="store")
requiredNamed.add_argument('--outfile', metavar="out_file", required=True, help="Output file containing assembly size and number of mapped and unmapped reads", action="store")
args = parser.parse_args()

flagstat = args.flagstat
scaffolds = args.scaffolds
flagstat = args.flagstat
outfile = args.outfile
is_tsv = True if flagstat[0].endswith(".tsv") else False

def get_read_count(sample, file_info):
    num_reads = 0
    file_path = os.path.join(os.getcwd(), "results", file_info[0], sample + file_info[1])
    with gzip.open(file_path, "rt") as fq_fh:
        for record in SeqIO.parse(fq_fh, "fastq"):
            num_reads += 1
    return num_reads


counts_dict = {}
mapped_dict = {}
for file in flagstat:
    if "_mapped_flagstat.txt" in file:
        sample = file.split("/")[-1].split("_mapped_flagstat.txt")[0]
    else:
        sys.exit(f"Trying to parse sample name out of {file} but it is not in the format expected <sample>_mapped_flagstat.txt.")
    with open(file, "r") as fh:
        if is_tsv:
            for line in fh:
                line_split = line.strip().split("\t")
                if "primary" in line_split[2]:
                    count = line_split[0]
                    mapped_dict[sample] = int(count)
        else:
            for line in fh:
                line_split = line.strip().split(" ")
                if line_split[3] == "primary":
                    count = line_split[0]
                    mapped_dict[sample] = int(count)  
        counts_dict[sample] = 2*int(get_read_count(sample, ["01_cleanreads","_R1_repaired.fastq.gz"]))

Sample_list = list()
MinContigLen_list = list()
MaxContigLen_list = list()
NumberOfContigs_list = list()
NumberOfContigsUn_list = list()
AssemblySize_list = list()
NumReads_list = list()
AssemblyMapped_list = list()
ProportionMapped_list = list()
ContigsN50_list = list()

for i in range(len(flagstat)):
    print(i)
    file = flagstat[i]
    Sample = file.split("/")[-1].split("_mapped_flagstat.txt")[0]
    print(Sample)
    Sample_list.append(Sample)
    
    ContigLensList = [int(n) for n in (os.popen("echo $(cat "+scaffolds[i]+" | grep \'>\' | cut -d\'_\' -f4 | tr \'\n\' \',\' )").read()).split(",") if n != '']
    tmp = []
    for tmp_number in set(ContigLensList):
            tmp += [tmp_number] * ContigLensList.count(tmp_number) * tmp_number
    tmp.sort()
    if (len(tmp) % 2) == 0:
        ContigsN50_list.append((tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2)
    else:
        ContigsN50_list.append(tmp[int(len(tmp) / 2)])

    MinContigLen = int(os.popen("echo $(cat "+scaffolds[i]+" | grep \'>\' | cut -d\'_\' -f4 | sort -n | head -1 )").read())
    MinContigLen_list.append(MinContigLen)

    MaxContigLen = int(os.popen("echo $(cat "+scaffolds[i]+" | grep \'>\' | cut -d\'_\' -f4 | sort -n | tail -1 )").read())
    MaxContigLen_list.append(MaxContigLen)

    NumberOfContigs = int(os.popen("echo $(cat "+scaffolds[i]+" | grep -c \'>\')").read())
    NumberOfContigs_list.append(NumberOfContigs)

    AssemblySize = int(os.popen("echo $(cat "+scaffolds[i]+" | grep -v \'>\' | tr -d \'\n\' | wc -m)").read())
    AssemblySize_list.append(AssemblySize)

    NumReads = counts_dict[Sample]
    NumReads_list.append(NumReads)

    AssemblyMapped = mapped_dict[Sample]
    AssemblyMapped_list.append(AssemblyMapped)

ProportionMapped_list = [str(round(x/y*100, 2)) for (x,y) in zip(AssemblyMapped_list, NumReads_list)]

Summary_path = os.path.join(os.getcwd(), outfile)

with open(Summary_path, "w") as file:
    file.write("Sample, Assembly size, Number of reads, " +
    "Number of filtered scaffolds, " +
    "N50, Min contig length, Max contig length, " +
    "Number mapped, Percent mapped\n")
    for i in range(len(input.reads1)):
        file.write(f"{Sample_list[i]}, {AssemblySize_list[i]}, {NumReads_list[i]}, {NumberOfContigs_list[i]}, {ContigsN50_list[i]}, {MinContigLen_list[i]}, {MaxContigLen_list[i]}, {AssemblyMapped_list[i]}, {ProportionMapped_list[i]}\n")