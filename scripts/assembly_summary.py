import os
import sys

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('--flagstat', nargs="+", metavar="flagstat",required=True, help="space-separated list of paths to flagstat files", action="store")
requiredNamed.add_argument('--scaffolds', nargs="+", metavar="scaffolds",required=True, help="space-separated list of paths to scaffolds files", action="store")
requiredNamed.add_argument('--scaffolds_unparsed', nargs="+", metavar="scaffolds_unparsed",required=True, help="space-separated list of paths to scaffolds_unparsed files", action="store")
requiredNamed.add_argument('--outfile', metavar="out_file", required=True, help="Output file containing assembly size and number of mapped and unmapped reads", action="store")
args = parser.parse_args()

flagstat = args.flagstat
scaffolds = args.scaffolds
scaffolds_unparsed = args.scaffolds_unparsed
flagstat = args.flagstat
outfile = args.outfile

counts_dict = {}
mapped_dict = {}
for file in flagstat:
    if "_assembly_mapping_flagstat.tsv" in file:
        sample = file.split("/")[-1].split("_assembly_mapping_flagstat.tsv")[0]
    else:
        sys.exit(f"Trying to parse sample name out of {file} but it is not in the format expected <sample>_assembly_mapping_flagstat.tsv.")
    with open(file, "r") as fh:
        for line in fh:
            line_split = line.strip().split("\t")
            if line_split[2] == "mapped":
                count = line_split[0]
                mapped_dict[sample] = int(count)
            if "primary" in line_split[2]:
                count = line_split[0]
                counts_dict[sample] = int(count)
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
    Sample = file.split("/")[-1].split("_assembly_mapping_flagstat.tsv")[0]
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

    NumberOfContigsUn = int(os.popen("echo $(cat "+scaffolds_unparsed[i]+" | grep -c \'>\')").read())
    NumberOfContigsUn_list.append(NumberOfContigsUn)

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
    "Total nummber of scaffolds, Number of filtered scaffolds," +
    "N50 filtered, Min contig length, Max contig length, " +
    "Number mapped, Percent mapped\n")
    for i in range(len(input.reads1)):
        file.write(f"{Sample_list[i]}, {AssemblySize_list[i]}, {NumReads_list[i]}, {NumberOfContigsUn_list[i]}, {NumberOfContigs_list[i]}, {ContigsN50_list[i]}, {MinContigLen_list[i]}, {MaxContigLen_list[i]}, {AssemblyMapped_list[i]}, {ProportionMapped_list[i]}\n")