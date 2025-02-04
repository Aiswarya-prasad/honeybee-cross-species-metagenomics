mkdir additional_analyses/results/add_isolates/g__Bombilactobacillus
mkdir additional_analyses/results/add_isolates/genomes_downloaded

# downloade assembly summary text
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O additional_analyses/results/add_isolates/assembly_summary.txt


# Download genomes

# Bombilactobacillus:
# GCF_000970795.1	PRJNA224116	SAMN03275731	JXJQ00000000.1	reference genome	1218492	1218492	Bombilactobacillus mellifer	strain=Bin4	na	latest	Scaffold	Major	Full	2015/04/13	ASM97079v1	Uppsala university	GCA_000970795.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/970/795/GCF_000970795.1_ASM97079v1	na	assembly from type material	na	haploid	bacteria	1815047	1786090	39.500000	11	28	NCBI RefSeq	GCF_000970795.1-RS_2024_05_13	2024-05-13	1734	1658	58	na
# GCF_025290075.1	PRJNA224116	SAMN27531748	JAMBKR000000000.1	na	1218492	1218492	Bombilactobacillus mellifer	na	MRS2-bin.14	latest	Scaffold	Major	Full	2022/09/19	ASM2529007v1	China Agricultural University	GCA_025290075.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/290/075/GCF_025290075.1_ASM2529007v1	derived from metagenome	na	na	haploid	bacteria	1655080	1655050	40.000000	0	67	67	NCBI RefSeq	GCF_025290075.1-RS_2024_03_29	2024-03-29	1638	1586	32	36045431
# GCF_025290985.1	PRJNA224116	SAMN27531741	JAMBKK000000000.1	na	1218492	1218492	Bombilactobacillus mellifer	na	MRS1-bin.8	latest	Scaffold	Major	Full	2022/09/19	ASM2529098v1	China Agricultural University	GCA_025290985.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/290/985/GCF_025290985.1_ASM2529098v1	derived from metagenome	na	na	haploid	bacteria	1523064	1523034	40.000000	0	48	48	NCBI RefSeq	GCF_025290985.1-RS_2024_03_29	2024-03-29	1500	1459	25	36045431
# GCF_025291455.1	PRJNA224116	SAMN27531715	JAMBJK000000000.1	na	1218492	1218492	Bombilactobacillus mellifer	na	GUT-bin.11	latest	Scaffold	Major	Full	2022/09/19	ASM2529145v1	China Agricultural University	GCA_025291455.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/291/455/GCF_025291455.1_ASM2529145v1	derived from metagenome	na	na	haploid	bacteria	1512989	1509570	40.000000	0	38	38	NCBI RefSeq	GCF_025291455.1-RS_2024_03_29	2024-03-29	1465	1441	15	36045431
# GCF_042663165.1	PRJNA224116	SAMN43283206	JBHSZT000000000.1	na	1218492	1218492	Bombilactobacillus mellifer	strain=CCUG 57507	na	latest	Contig	Major	Full	2024/10/06	ASM4266316v1	Institution of Microbiology, Chinese Academy of Science(WDCM)	GCA_042663165.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/663/165/GCF_042663165.1_ASM4266316v1	na	na	na	haploid	bacteria	1790866	1790866	39.500000	0	10	10	NCBI RefSeq	GCF_042663165.1-RS_2024_10_06	2024-10-06	1742	1663	65	30832757
# GCF_000967245.1	PRJNA224116	SAMN03271966	JXBZ00000000.1	reference genome	1218508	1218508	Bombilactobacillus mellis	strain=Hon2	na	latest	Scaffold	Major	Full	2015/03/31	ASM96724v1	Uppsala university	GCA_000967245.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/967/245/GCF_000967245.1_ASM96724v1	na	assembly from type material	na	haploid	bacteria	1810599	1790038	36.000000	17	NCBI RefSeq	GCF_000967245.1-RS_2024_05_13	2024-05-13	1724	1652	59	na
# GCF_013345055.1	PRJNA224116	SAMN13893454	JAAEDZ000000000.1	na	1218508	1218508	Bombilactobacillus mellis	strain=ESL0295	na	latest	Contig	Major	Full	2020/06/15	ASM1334505v1	University of LausannGCA_013345055.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/345/055/GCF_013345055.1_ASM1334505v1	na	na	na	haploid	bacteria	1683724	1683724	36.500000	0	12	12	NCBI RefSeq	GCF_013345055.1-RS_2024_03_29	2024-03-29	1627	1560	55	32531278
# GCF_013346925.1	PRJNA224116	SAMN13893453	JAAEEA000000000.1	na	1218508	1218508	Bombilactobacillus mellis	strain=ESL0294	na	latest	Contig	Major	Full	2020/06/15	ASM1334692v1	University of LausannGCA_013346925.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/346/925/GCF_013346925.1_ASM1334692v1	na	na	na	haploid	bacteria	1729063	1729063	36.500000	0	17	17	NCBI RefSeq	GCF_013346925.1-RS_2024_03_29	2024-03-29	1644	1570	63	32531278
# GCF_013347085.1	PRJNA224116	SAMN13893455	JAAEDY000000000.1	na	1218508	1218508	Bombilactobacillus mellis	strain=ESL0394	na	latest	Contig	Major	Full	2020/06/15	ASM1334708v1	University of LausannGCA_013347085.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/347/085/GCF_013347085.1_ASM1334708v1	na	na	na	haploid	bacteria	1743889	1743889	36.500000	0	16	16	NCBI RefSeq	GCF_013347085.1-RS_2024_03_29	2024-03-29	1664	1577	69	32531278
# GCF_026229285.1	PRJNA224116	SAMN30153779	JANUHE000000000.1	na	1218508	1218508	Bombilactobacillus mellis	strain=LB26	na	latest	Contig	Major	Full	2022/11/17	ASM2622928v1	Virginia Tech	GCA_026229285.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/229/285/GCF_026229285.1_ASM2622928v1	na	na	na	haploid	bacteria	1644415	1644415	36.500000	0	15	15	NCBI RefSeq	NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	2024-01-24	1550	1506	40	36331337
# GCF_042676405.1	PRJNA224116	SAMN43284444	JBHTHW000000000.1	na	1218508	1218508	Bombilactobacillus mellis	strain=CCUG 63289	na	latest	Contig	Major	Full	2024/10/06	ASM4267640v1	Institution of Microbiology, Chinese Academy of Science(WDCM)	GCA_042676405.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/676/405/GCF_042676405.1_ASM4267640v1	na	assembly from type material	na	haploid	bacteria	1789619	1789619	36.000000	0	11	11	NCBI RefSeq	GCF_042676405.1-RS_2024_10_07	2024-10-07	1719	1649	54	30832757
# GCF_003515755.1	PRJNA224116	SAMN09629877	QOCS00000000.1	na	1303590	1303590	Bombilactobacillus bombi	strain=LV-8.1	na	latest	Scaffold	Major	Full	2018/09/10	ASM351575v1	University of Texas at Austin	GCA_003515755.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/515/755/GCF_003515755.1_ASM351575v1	na	na	na	haploid	bacteria	1934723	1934200	34.500000	0	42	42	NCBI RefSeq	GCF_003515755.1-RS_2024_09_18	2024-09-18	2013	1876	59	na
# GCF_003515805.1	PRJNA224116	SAMN09629876	QOCR00000000.1	na	1303590	1303590	Bombilactobacillus bombi	strain=BI-1.1	na	latest	Contig	Major	Full	2018/09/10	ASM351580v1	University of Texas at AustinGCA_003515805.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/515/805/GCF_003515805.1_ASM351580v1	na	na	na	haploid	bacteria	1903140	1903140	34.500000	0	9	9	NCBI RefSeq	GCF_003515805.1-RS_2024_09_18	2024-09-18	1872	1801	57	na
# GCF_003522965.1	PRJNA224116	SAMN09629875	na	reference genome	1303590	1303590	Bombilactobacillus bombi	strain=BI-2.5	na	latest	Complete Genome	Major	Full	2018/09/10	ASM352296v1	University of Texas at Austin	GCA_003522965.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/522/965/GCF_003522965.1_ASM352296v1	na	na	na	haploid	bacteria	1842084	1842084	34.500000	1	1	1	NCBI RefSeq	GCF_003522965.1-RS_2024_03_28	2024-03-28	1768	1679	73	na
# GCF_013607485.1	PRJNA224116	SAMN10754285	SCHP00000000.1	na	1303590	1303590	Bombilactobacillus bombi	strain=XV6	na	latest	Contig	Major	Full	2020/07/23	ASM1360748v1	University of Bologna	GCA_013607485.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/607/485/GCF_013607485.1_ASM1360748v1	na	na	na	haploid	bacteria	1805736	1805736	36.000000	0	164	164	NCBI RefSeq	GCF_013607485.1-RS_2024_09_07	2024-09-07	1811	1715	49	na
# GCF_042432605.1	PRJNA224116	SAMN43281635	JBHLWZ000000000.1	na	1303590	1303590	Bombilactobacillus bombi	strain=CCM 8440	na	latest	Contig	Major	Full	2024/09/28	ASM4243260v1	Institution of Microbiology, Chinese Academy of Science(WDCM)	GCA_042432605.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/432/605/GCF_042432605.1_ASM4243260v1	na	assembly from type material	na	haploid	bacteria	1909329	1909329	34.500000	0	22	22	NCBI RefSeq	GCF_042432605.1-RS_2024_09_29	2024-09-29	1835	1753	62	na
# GCF_013385145.1	PRJNA224116	SAMN15312313	JABZEC000000000.1	reference genome	2675299	2675299	Bombilactobacillus apium	strain=DCY120	na	latest	Scaffold	Major	Full	2020/07/01	ASM1338514v1	Kyung Hee University	GCA_013385145.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/385/145/GCF_013385145.1_ASM1338514v1	na	assembly from type material	na	haploid	bacteria	1712546	1712451	40.000000	0	16	16	NCBI RefSeq	NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	2024-01-11	1747	1630	62	na
# GCF_023380265.1	PRJNA224116	SAMN26370965	na	reference genome	2923362	2923362	Bombilactobacillus folatiphilus	strain=SG4_D2	na	latest	Complete Genome	Major	Full	2022/05/17	ASM2338026v1	The University of Adelaide	GCA_023380265.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/380/265/GCF_023380265.1_ASM2338026v1	na	assembly from type material	na	haploid	bacteria	1637944	1637944	38.500000	NCBI RefSeq	NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	2023-10-29	1627	1536	73	36094463
# GCF_023380245.1	PRJNA224116	SAMN26370966	na	reference genome	2923363	2923363	Bombilactobacillus thymidiniphilus	strain=SG4_A1	na	latest	Complete Genome	Major	Full	2022/05/17	ASM2338024v1	The University of Adelaide	GCA_023380245.1	identical	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/380/245/GCF_023380245.1_ASM2338024v1	na	assembly from type material	na	haploid	bacteria	1494436	1494436	36.500000	1	1	1	NCBI RefSeq	NCBI Prokaryotic Genome Annotation Pipeline (PGAP)	2023-10-29	1535	1449	73	36094463

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/970/795/GCF_000970795.1_ASM97079v1/GCF_000970795.1_ASM97079v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_000970795.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/290/075/GCF_025290075.1_ASM2529007v1/GCF_025290075.1_ASM2529007v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_025290075.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/290/985/GCF_025290985.1_ASM2529098v1/GCF_025290985.1_ASM2529098v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_025290985.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/291/455/GCF_025291455.1_ASM2529145v1/GCF_025291455.1_ASM2529145v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_025291455.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/663/165/GCF_042663165.1_ASM4266316v1/GCF_042663165.1_ASM4266316v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_042663165.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/967/245/GCF_000967245.1_ASM96724v1/GCF_000967245.1_ASM96724v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_000967245.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/345/055/GCF_013345055.1_ASM1334505v1/GCF_013345055.1_ASM1334505v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_013345055.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/346/925/GCF_013346925.1_ASM1334692v1/GCF_013346925.1_ASM1334692v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_013346925.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/347/085/GCF_013347085.1_ASM1334708v1/GCF_013347085.1_ASM1334708v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_013347085.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/229/285/GCF_026229285.1_ASM2622928v1/GCF_026229285.1_ASM2622928v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_026229285.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/676/405/GCF_042676405.1_ASM4267640v1/GCF_042676405.1_ASM4267640v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_042676405.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/515/755/GCF_003515755.1_ASM351575v1/GCF_003515755.1_ASM351575v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_003515755.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/515/805/GCF_003515805.1_ASM351580v1/GCF_003515805.1_ASM351580v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_003515805.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/522/965/GCF_003522965.1_ASM352296v1/GCF_003522965.1_ASM352296v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_003522965.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/607/485/GCF_013607485.1_ASM1360748v1/GCF_013607485.1_ASM1360748v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_013607485.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/432/605/GCF_042432605.1_ASM4243260v1/GCF_042432605.1_ASM4243260v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_042432605.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/385/145/GCF_013385145.1_ASM1338514v1/GCF_013385145.1_ASM1338514v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_013385145.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/380/265/GCF_023380265.1_ASM2338026v1/GCF_023380265.1_ASM2338026v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_023380265.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/380/245/GCF_023380245.1_ASM2338024v1/GCF_023380245.1_ASM2338024v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/GCF_023380245.1.fna.gz


gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_000970795.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_025290075.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_025290985.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_025291455.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_042663165.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_000967245.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_013345055.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_013346925.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_013347085.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_026229285.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_042676405.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_003515755.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_003515805.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_003522965.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_013607485.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_042432605.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_013385145.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_023380265.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/GCF_023380245.1.fna.gz

# annotate with prodigal
conda activate /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs/c41bb15eb344238b9c5a180c961ff099_
# PRODIGAL v2.6.3 [February, 2016]         
# Univ of Tenn / Oak Ridge National Lab
# Doug Hyatt, Loren Hauser, et al.     
# -------------------------------------
for genome in additional_analyses/results/add_isolates/genomes_downloaded/*.fna
do
    genome_name=$(basename ${genome} .fna)
    genome_gff=additional_analyses/results/add_isolates/genomes_downloaded/${genome_name}.gff
    prodigal -i ${genome} -a additional_analyses/results/add_isolates/genomes_downloaded/${genome_name}.faa -o ${genome_gff} -f gff
done

# rename such that each header is genomeaccession_n where n is a sequential number starting from 1
for genome in additional_analyses/results/add_isolates/genomes_downloaded/*.faa
do
    # rename the file as _original_headers
    mv ${genome} ${genome/.faa/_original_headers.faa}
    # mv ${genome} ${genome/.gff/_original_headers.gff}
    # mv ${genome} ${genome/.fna/_original_headers.fna}
done

for genome in additional_analyses/results/add_isolates/genomes_downloaded/*_original_headers.faa
do
    genome_name=$(basename ${genome})
    genome_name=${genome_name/_original_headers.faa/}
    # add the genome_name before serial number
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/genomes_downloaded/${genome_name}.faa
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/g__Bombilactobacillus/isolates_to_add/${genome_name}.faa
done

# copy orthofinder working dir
rsync -aviP /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/11_phylogenies/02_orthofinder_results/g__Bombilactobacillus/Results_g__Bombilactobacillus_iqtree additional_analyses/results/add_isolates/g__Bombilactobacillus/Results_g__Bombilactobacillus_iqtree

sbatch additional_analyses/scripts/run_orthofinder_job_g__Bombilactobacillus.sh


# download genomes for other Lactobacillus into additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus
mkdir -p additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/837/055/GCF_002837055.1_ASM283705v1/GCF_002837055.1_ASM283705v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_002837055.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/112/665/GCF_900112665.1_IMG-taxon_2617270727_annotated_assembly/GCF_900112665.1_IMG-taxon_2617270727_annotated_assembly_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_900112665.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/428/255/GCF_026428255.1_ASM2642825v1/GCF_026428255.1_ASM2642825v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_026428255.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/972/815/GCF_019972815.1_ASM1997281v1/GCF_019972815.1_ASM1997281v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_019972815.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/970/755/GCF_000970755.1_ASM97075v1/GCF_000970755.1_ASM97075v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_000970755.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/967/195/GCF_000967195.1_ASM96719v1/GCF_000967195.1_ASM96719v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_000967195.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/972/835/GCF_019972835.1_ASM1997283v1/GCF_019972835.1_ASM1997283v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_019972835.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/185/255/GCF_026185255.2_ASM2618525v2/GCF_026185255.2_ASM2618525v2_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_026185255.2.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/693/045/GCF_003693045.1_ASM369304v1/GCF_003693045.1_ASM369304v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_003693045.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/469/265/GCF_019469265.1_ASM1946926v1/GCF_019469265.1_ASM1946926v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_019469265.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/100/975/GCF_016100975.1_ASM1610097v1/GCF_016100975.1_ASM1610097v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_016100975.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/323/605/GCF_014323605.1_ASM1432360v1/GCF_014323605.1_ASM1432360v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_014323605.1.fna.gz

gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_002837055.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_900112665.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_026428255.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_019972815.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_000970755.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_000967195.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_019972835.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_026185255.2.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_003693045.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_019469265.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_016100975.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/GCF_014323605.1.fna.gz

# annotate with prodigal
conda activate /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs/c41bb15eb344238b9c5a180c961ff099_
# PRODIGAL v2.6.3 [February, 2016]
# Univ of Tenn / Oak Ridge National Lab
# Doug Hyatt, Loren Hauser, et al.
# -------------------------------------
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/*.fna
do
    genome_name=$(basename ${genome} .fna)
    genome_gff=additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/${genome_name}.gff
    prodigal -i ${genome} -a additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/${genome_name}.faa -o ${genome_gff} -f gff
done

# rename such that each header is genomeaccession_n where n is a sequential number starting from 1
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/*.faa
do
    # rename the file as _original_headers
    mv ${genome} ${genome/.faa/_original_headers.faa}
    # mv ${genome} ${genome/.gff/_original_headers.gff}
    # mv ${genome} ${genome/.fna/_original_headers.fna}
done

mkdir -p additional_analyses/results/add_isolates/g__Lactobacillus/
mkdir -p additional_analyses/results/add_isolates/g__Lactobacillus/isolates_to_add/
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/*_original_headers.faa
do
    genome_name=$(basename ${genome})
    genome_name=${genome_name/_original_headers.faa/}
    # add the genome_name before serial number
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/genomes_downloaded/g__Lactobacillus/${genome_name}.faa
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/g__Lactobacillus/isolates_to_add/${genome_name}.faa
done

# copy orthofinder working dir
rsync -aviP /work/FAC/FBM/DMF/pengel/spirit/aprasad/BACKUP_current/20230313_apis_species_comparison/results/11_phylogenies/02_orthofinder_results/g__Lactobacillus/Results_g__Lactobacillus additional_analyses/results/add_isolates/g__Lactobacillus/Results_g__Lactobacillus

sbatch additional_analyses/scripts/run_orthofinder_job_g__Lactobacillus.sh



# download genomes for other Bifidobacterium into additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium
mkdir -p additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/007/559/275/GCF_007559275.1_ASM755927v1/GCF_007559275.1_ASM755927v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_007559275.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/715/865/GCF_002715865.1_ASM271586v1/GCF_002715865.1_ASM271586v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_002715865.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/202/755/GCF_003202755.1_ASM320275v1/GCF_003202755.1_ASM320275v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_003202755.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/202/695/GCF_003202695.1_ASM320269v1/GCF_003202695.1_ASM320269v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_003202695.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/683/175/GCF_009683175.1_ASM968317v1/GCF_009683175.1_ASM968317v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_009683175.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/469/425/GCF_019469425.1_ASM1946942v1/GCF_019469425.1_ASM1946942v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_019469425.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/102/005/GCF_016102005.1_ASM1610200v1/GCF_016102005.1_ASM1610200v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_016102005.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/706/765/GCF_000706765.1_ASM70676v1/GCF_000706765.1_ASM70676v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_000706765.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/884/755/GCF_020884755.1_ASM2088475v1/GCF_020884755.1_ASM2088475v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_020884755.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/101/585/GCF_016101585.1_ASM1610158v1/GCF_016101585.1_ASM1610158v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_016101585.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/499/285/GCF_000499285.1_bifidobacterium_strain7101/GCF_000499285.1_bifidobacterium_strain7101_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_000499285.1.fna.gz

gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_007559275.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_002715865.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_003202755.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_003202695.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_009683175.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_019469425.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_016102005.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_000706765.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_020884755.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_016101585.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/GCF_000499285.1.fna.gz

# annotate with prodigal
conda activate /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs/c41bb15eb344238b9c5a180c961ff099_
# PRODIGAL v2.6.3 [February, 2016]
# Univ of Tenn / Oak Ridge National Lab
# Doug Hyatt, Loren Hauser, et al.
# -------------------------------------
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/*.fna
do
    genome_name=$(basename ${genome} .fna)
    genome_gff=additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/${genome_name}.gff
    prodigal -i ${genome} -a additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/${genome_name}.faa -o ${genome_gff} -f gff
done

# rename such that each header is genomeaccession_n where n is a sequential number starting from 1
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/*.faa
do
    # rename the file as _original_headers
    mv ${genome} ${genome/.faa/_original_headers.faa}
    # mv ${genome} ${genome/.gff/_original_headers.gff}
    # mv ${genome} ${genome/.fna/_original_headers.fna}
done

mkdir -p additional_analyses/results/add_isolates/g__Bifidobacterium/
mkdir -p additional_analyses/results/add_isolates/g__Bifidobacterium/isolates_to_add/
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/*_original_headers.faa
do
    genome_name=$(basename ${genome})
    genome_name=${genome_name/_original_headers.faa/}
    # add the genome_name before serial number
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/genomes_downloaded/g__Bifidobacterium/${genome_name}.faa
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/g__Bifidobacterium/isolates_to_add/${genome_name}.faa
done

# copy orthofinder working dir
rsync -aviP results/11_phylogenies/02_orthofinder_results/g__Bifidobacterium/Results_g__Bifidobacterium additional_analyses/results/add_isolates/g__Bifidobacterium

sbatch additional_analyses/scripts/run_orthofinder_job_g__Bifidobacterium.sh


# download genomes for other Gilliamella into additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella
mkdir -p additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/599/985/GCF_000599985.1_ASM59998v1/GCF_000599985.1_ASM59998v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_000599985.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/469/165/GCF_019469165.1_ASM1946916v1/GCF_019469165.1_ASM1946916v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_019469165.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/693/755/GCF_001693755.1_ASM169375v1/GCF_001693755.1_ASM169375v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001693755.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/690/185/GCF_001690185.1_ASM169018v1/GCF_001690185.1_ASM169018v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001690185.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/690/685/GCF_001690685.1_ASM169068v1/GCF_001690685.1_ASM169068v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001690685.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/690/705/GCF_001690705.1_ASM169070v1/GCF_001690705.1_ASM169070v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001690705.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/690/195/GCF_001690195.1_ASM169019v1/GCF_001690195.1_ASM169019v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001690195.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/693/435/GCF_001693435.1_ASM169343v1/GCF_001693435.1_ASM169343v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001693435.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/202/915/GCF_003202915.1_ASM320291v1/GCF_003202915.1_ASM320291v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_003202915.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/142/155/GCF_002142155.1_ASM214215v1/GCF_002142155.1_ASM214215v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_002142155.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/758/615/GCF_030758615.1_ASM3075861v1/GCF_030758615.1_ASM3075861v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_030758615.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/142/215/GCF_002142215.1_ASM214221v1/GCF_002142215.1_ASM214221v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_002142215.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/469/185/GCF_019469185.1_ASM1946918v1/GCF_019469185.1_ASM1946918v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_019469185.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/469/205/GCF_019469205.1_ASM1946920v1/GCF_019469205.1_ASM1946920v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_019469205.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/751/545/GCF_028751545.1_ASM2875154v1/GCF_028751545.1_ASM2875154v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_028751545.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/536/085/GCF_026536085.1_ASM2653608v1/GCF_026536085.1_ASM2653608v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_026536085.1.fna.gz

gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_000599985.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_019469165.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001693755.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001690185.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001690685.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001690705.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001690195.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_001693435.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_003202915.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_002142155.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_030758615.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_002142215.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_019469185.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_019469205.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_028751545.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/GCF_026536085.1.fna.gz

# annotate with prodigal
conda activate /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs/c41bb15eb344238b9c5a180c961ff099_
# PRODIGAL v2.6.3 [February, 2016]
# Univ of Tenn / Oak Ridge National Lab
# Doug Hyatt, Loren Hauser, et al.
# -------------------------------------
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/*.fna
do
    genome_name=$(basename ${genome} .fna)
    genome_gff=additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/${genome_name}.gff
    prodigal -i ${genome} -a additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/${genome_name}.faa -o ${genome_gff} -f gff
done

# rename such that each header is genomeaccession_n where n is a sequential number starting from 1
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/*.faa
do
    # rename the file as _original_headers
    mv ${genome} ${genome/.faa/_original_headers.faa}
    # mv ${genome} ${genome/.gff/_original_headers.gff}
    # mv ${genome} ${genome/.fna/_original_headers.fna}
done

mkdir -p additional_analyses/results/add_isolates/g__Gilliamella/
mkdir -p additional_analyses/results/add_isolates/g__Gilliamella/isolates_to_add/
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/*_original_headers.faa
do
    genome_name=$(basename ${genome})
    genome_name=${genome_name/_original_headers.faa/}
    # add the genome_name before serial number
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/genomes_downloaded/g__Gilliamella/${genome_name}.faa
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/g__Gilliamella/isolates_to_add/${genome_name}.faa
done

# copy orthofinder working dir
rsync -aviP results/11_phylogenies/02_orthofinder_results/g__Gilliamella/Results_g__Gilliamella additional_analyses/results/add_isolates/g__Gilliamella

sbatch additional_analyses/scripts/run_orthofinder_job_g__Gilliamella.sh



# download genomes for other Snodgrassella into additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella

mkdir -p additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/600/005/GCF_000600005.1_ASM60000v1/GCF_000600005.1_ASM60000v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/GCF_000600005.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/777/855/GCF_002777855.1_ASM277785v1/GCF_002777855.1_ASM277785v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/GCF_002777855.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/777/745/GCF_002777745.1_ASM277774v1/GCF_002777745.1_ASM277774v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/GCF_002777745.1.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/535/915/GCF_026535915.1_ASM2653591v1/GCF_026535915.1_ASM2653591v1_genomic.fna.gz -O additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/GCF_026535915.1.fna.gz

gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/GCF_000600005.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/GCF_002777855.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/GCF_002777745.1.fna.gz
gunzip additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/GCF_026535915.1.fna.gz

# annotate with prodigal
conda activate /work/FAC/FBM/DMF/pengel/spirit/aprasad/snakemake-conda-envs/c41bb15eb344238b9c5a180c961ff099_
# PRODIGAL v2.6.3 [February, 2016]
# Univ of Tenn / Oak Ridge National Lab
# Doug Hyatt, Loren Hauser, et al.
# -------------------------------------
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/*.fna
do
    genome_name=$(basename ${genome} .fna)
    genome_gff=additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/${genome_name}.gff
    prodigal -i ${genome} -a additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/${genome_name}.faa -o ${genome_gff} -f gff
done

# rename such that each header is genomeaccession_n where n is a sequential number starting from 1
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/*.faa
do
    # rename the file as _original_headers
    mv ${genome} ${genome/.faa/_original_headers.faa}
    # mv ${genome} ${genome/.gff/_original_headers.gff}
    # mv ${genome} ${genome/.fna/_original_headers.fna}
done

mkdir -p additional_analyses/results/add_isolates/g__Snodgrassella/
mkdir -p additional_analyses/results/add_isolates/g__Snodgrassella/isolates_to_add/
for genome in additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/*_original_headers.faa
do
    genome_name=$(basename ${genome})
    genome_name=${genome_name/_original_headers.faa/}
    # add the genome_name before serial number
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/genomes_downloaded/g__Snodgrassella/${genome_name}.faa
    awk -v name="${genome_name}" '/^>/{print ">" name "_" ++i; next}{print}' < ${genome} > additional_analyses/results/add_isolates/g__Snodgrassella/isolates_to_add/${genome_name}.faa
done

# copy orthofinder working dir
rsync -aviP results/11_phylogenies/02_orthofinder_results/g__Snodgrassella/Results_g__Snodgrassella additional_analyses/results/add_isolates/g__Snodgrassella

sbatch additional_analyses/scripts/run_orthofinder_job_g__Snodgrassella.sh