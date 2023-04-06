#!/bin/bash

prefix_source="aiswarya@130.223.110.124:/home/aiswarya/mnt/nas/ToBeTransferred_nas_recherche/datasets"
prefix_destination="/scratch/aprasad/20230313_apis_species_comparison/results"

declare -A names=(
[Ig13612]=AmAi02
[Ig13618]=AmAi08
[Ig13631]=AcCh01
[Ig13632]=AcCh02
[Ig13611]=AmAi01
[Ig13613]=AmAi03
[Ig13614]=AmAi04
[Ig13615]=AmAi05
[Ig13616]=AmAi06
[Ig13617]=AmAi07
[Ig13619]=AmAi09
[Ig13620]=AmAi10
[Ig13621]=AmIu01
[Ig13622]=AmIu02
[Ig13623]=AmIu03
[Ig13624]=AmIu04
[Ig13625]=AmIu05
[Ig13626]=AmIu06
[Ig13627]=AmIu07
[Ig13628]=AmIu08
[Ig13629]=AmIu09
[Ig13630]=AmIu10
[Ig13633]=AcCh03
[Ig13634]=AcCh04
[Ig13635]=AcCh05
[Ig13636]=AcCh06
[Ig13637]=AcCh07
[Ig13638]=AcCh08
[Ig13639]=AcCh09
[Ig13640]=AcCh10
[Ig13641]=AcKn01
[Ig13642]=AcKn02
[Ig13643]=AcKn03
[Ig13644]=AcKn04
[Ig13645]=AcKn05
[Ig13646]=AcKn06
[Ig13647]=AcKn07
[Ig13648]=AcKn08
[Ig13649]=AcKn09
[Ig13650]=AcKn10
);

for name in "${!names[@]}"
do
  codename="$name"
  finalname="${names[$name]}"
  echo  "copying $codename and calling it $finalname"
  rsync -avi --progress ${prefix_source}/NGS_data/20180612_KE_japan_metagenomes/data/"$codename"_R1.fastq.gz ${prefix_destination}/NGS_data/20180612_KE_japan_metagenomes/
  rsync -avi --progress ${prefix_source}/NGS_data/20180612_KE_japan_metagenomes/data/"$codename"_R2.fastq.gz ${prefix_destination}/NGS_data/20180612_KE_japan_metagenomes/
done


rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F1_R1.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F1_R2.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F2_R1.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F2_R2.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F3_R1.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F3_R2.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F4_R1.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F4_R2.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F5_R1.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/
rsync -avi --progress ${prefix_source}/NGS_data/20170310_WINDU179/data/DrY2_F5_R2.fastq.gz ${prefix_destination}/NGS_data/20170310_WINDU179/


rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N1_R1.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N1_R2.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N2_R1.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N2_R2.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N3_R1.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N3_R2.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N4_R1.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N4_R2.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N5_R1.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N5_R2.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N6_R1.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/
rsync -avi --progress ${prefix_source}/NGS_data/20151119_WINDU89/data/DrY1_N6_R2.fastq.gz ${prefix_destination}/NGS_data/20151119_WINDU89/


rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F1_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F1_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F2_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F2_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F3_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F3_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F4_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F4_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F5_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F5_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F6_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_F6_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W1_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W1_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W2_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W2_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W3_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W3_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W4_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W4_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W5_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W5_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W6_R1.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/
rsync -avi --progress ${prefix_source}/NGS_data/20160415_OBIWAN225/DrY1_W6_R2.fastq.gz ${prefix_destination}/NGS_data/20160415_OBIWAN225/


rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N1_R1.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N1_R2.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N2_R1.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N2_R2.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N3_R1.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N3_R2.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N4_R1.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N4_R2.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N5_R1.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N5_R2.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N6_R1.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/
rsync -avi --progress ${prefix_source}/NGS_data/20161216_OBIWAN275/DrY2_N6_R2.fastq.gz ${prefix_destination}/NGS_data/20161216_OBIWAN275/


rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N1_R1.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N1_R2.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N2_R1.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N2_R2.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N3_R1.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N3_R2.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N4_R1.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N4_R2.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N5_R1.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N5_R2.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N6_R1.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/
rsync -avi --progress ${prefix_source}/NGS_data/20170426_OBIWAN300/GrY2_N6_R2.fastq.gz ${prefix_destination}/NGS_data/20170426_OBIWAN300/


rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F1_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F1_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F2_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F2_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F3_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F3_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F4_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F4_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F5_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F5_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F6_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_F6_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W1_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W1_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W2_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W2_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W3_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W3_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W4_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W4_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W5_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W5_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W6_R1.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 
rsync -avi --progress ${prefix_source}/NGS_data/20170428_WINDU191/GrY2_W6_R2.fastq.gz ${prefix_destination}/NGS_data/20170428_WINDU191/ 