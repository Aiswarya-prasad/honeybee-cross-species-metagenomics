# at the moment also using as a notebook for info
########### antiSMASH 6.1.1 #############
# on KE workstation
# conda activate antismash



# script to be copied and run on the lab workstation where antismash is set up within conda and working
# divide into four parts so others can run simultaneously
# A1-1 C5-1 D6-1 F7-1 A1-2 C5-2 D6-2 F7-2 A1-3 C5-3 D6-3 F7-3 A1-4 C5-4 D6-4 F7-4 A1-5 C5-5 D6-5 F7-5 A2-1 C6-1 D7-1 F8-1 A2-2 C6-2 D7-2 F8-2 A2-3 C6-3 D7-3 F8-3 A2-4 C6-4 D7-4 F8-4 A2-5 C6-5 D7-5 F8-5 A3-1 C7-1 D8-1 F9-1 A3-2 C7-2 D8-2 F9-2 A3-3 C7-3 D8-3 F9-3 A3-4
# C7-4 D8-4 F9-4 A3-5 C7-5 D8-5 F9-5 A4-1 C8-1 D9-1 M1-1 A4-2 C8-2 D9-2 M1-2 A4-3 C8-3 D9-3 M1-3 A4-4 C8-4 D9-4 M1-4 A4-5 C8-5 D9-5 M1-5 A5-1 C9-1 F1-1 M2-1 A5-2 C9-2 F1-2 M2-2 A5-3 C9-3 F1-3 M2-3 A5-4 C9-4 F1-4 M2-4 A5-5 C9-5 F1-5 M2-5
# A6-1 D1-1 F2-1 M3-1 A6-2 D1-2 F2-2 M3-2 A6-3 D1-3 F2-3 M3-3 A6-4 D1-4 F2-4 M3-4 A6-5 D1-5 F2-5 M3-5 C1-1 D2-1 F3-1 M4-1 C1-2 D2-2 F3-2 M4-2 C1-3 D2-3 F3-3 M4-3 C1-4 D2-4 F3-4 M4-4 C1-5 D2-5 F3-5 M4-5 C2-1 D3-1 F4-1 M5-1 C2-2 D3-2 F4-2
# M5-2 C2-3 D3-3 F4-3 M5-3 C2-4 D3-4 F4-4 M5-4 C2-5 D3-5 F4-5 M5-5 C3-1 D4-1 F5-1 M6-1 C3-2 D4-2 F5-2 M6-2 C3-3 D4-3 F5-3 M6-3 C3-4 D4-4 F5-4 M6-4 C3-5 D4-5 F5-5 M6-5 C4-1 D5-1 F6-1 M7-1 C4-2 D5-2 F6-2 M7-2 C4-3 D5-3 F6-3 M7-3 C4-4 D5-4 F6-4 M7-4 C4-5 D5-5 F6-5 M7-5
cd /home/aiswarya/20230313_apis_species_comparison
conda activate antismash
# # screen antismash
# for sample in A1-1;
# screen antismash1
for sample in C5-1 D6-1 F7-1 A1-2 C5-2 D6-2 F7-2 A1-3 C5-3 D6-3 F7-3 A1-4 C5-4 D6-4 F7-4 A1-5 C5-5 D6-5 F7-5 A2-1 C6-1 D7-1 F8-1 A2-2 C6-2 D7-2 F8-2 A2-3 C6-3 D7-3 F8-3 A2-4 C6-4 D7-4 F8-4 A2-5 C6-5 D7-5 F8-5 A3-1 C7-1 D8-1 F9-1 A3-2 C7-2 D8-2 F9-2 A3-3 C7-3 D8-3 F9-3 A3-4;
# # screen antismash2 - running
# for sample in C7-4 D8-4 F9-4 A3-5 C7-5 D8-5 F9-5 A4-1 C8-1 D9-1 M1-1 A4-2 C8-2 D9-2 M1-2 A4-3 C8-3 D9-3 M1-3 A4-4 C8-4 D9-4 M1-4 A4-5 C8-5 D9-5 M1-5 A5-1 C9-1 F1-1 M2-1 A5-2 C9-2 F1-2 M2-2 A5-3 C9-3 F1-3 M2-3 A5-4 C9-4 F1-4 M2-4 A5-5 C9-5 F1-5 M2-5;
# # screen antismash3 - running
# for sample in A6-1 D1-1 F2-1 M3-1 A6-2 D1-2 F2-2 M3-2 A6-3 D1-3 F2-3 M3-3 A6-4 D1-4 F2-4 M3-4 A6-5 D1-5 F2-5 M3-5 C1-1 D2-1 F3-1 M4-1 C1-2 D2-2 F3-2 M4-2 C1-3 D2-3 F3-3 M4-3 C1-4 D2-4 F3-4 M4-4 C1-5 D2-5 F3-5 M4-5 C2-1 D3-1 F4-1 M5-1 C2-2 D3-2 F4-2;
# # screen antismash4 - running
# for sample in M5-2 C2-3 D3-3 F4-3 M5-3 C2-4 D3-4 F4-4 M5-4 C2-5 D3-5 F4-5 M5-5 C3-1 D4-1 F5-1 M6-1 C3-2 D4-2 F5-2 M6-2 C3-3 D4-3 F5-3 M6-3 C3-4 D4-4 F5-4 M6-4 C3-5 D4-5 F5-5 M6-5 C4-1 D5-1 F6-1 M7-1 C4-2 D5-2 F6-2 M7-2 C4-3 D5-3 F6-3 M7-3 C4-4 D5-4 F6-4 M7-4 C4-5 D5-5 F6-5 M7-5;
do
if [ -f results/15_antismash/${sample}/${sample}.gbk ]; then
  echo skipping sample $sample
  continue
fi
cat results/08_gene_content/01_profiling_bowtie2/${sample}_filt_genes.gff | \
 sed -e "s/NODE/${sample}_NODE/" > results/08_gene_content/01_profiling_bowtie2/${sample}_filt_genes_renamed_temp_for_antismash.gff
echo "removing previous results"
rm -rf results/15_antismash/${sample}
echo running sample $sample
antismash results/07_MAG_binng_QC/00_assembled_scaffolds/${sample}/${sample}_scaffolds.fasta \
  -c 8 --taxon bacteria --fullhmmer --clusterhmmer \
  --asf --cc-mibig --cb-general --cb-subclusters \
  --cb-knownclusters --pfam2go --rre \
  --output-dir results/15_antismash/${sample} \
  --output-basename ${sample} \
  --genefinding-gff3 results/08_gene_content/01_profiling_bowtie2/${sample}_filt_genes_renamed_temp_for_antismash.gff &> results/15_antismash/${sample}_antismash.log
done

# rsync -aviP 15_antismash/ aprasad@curnagl.dcsr.unil.ch:/work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/08_gene_content