prefix="aprasad@curnagl.dcsr.unil.ch:/scratch/aprasad/211018_Medgenome_india_samples"
outdir="/Volumes/Storage/Work/Temp-From-NAS/cross-species-analysis"
# checkpoints
mkdir -p 01_Trimmed/
rsync -r --progress -v $prefix/01_Trimmed/*_R*_trim.fastq.gz $outdir/01_Trimmed/*_R*_trim.fastq.gz
mkdir -p 05_Assembly/
mkdir -p 05_Assembly/host_unmapped
rsync -r --progress -v $prefix/05_Assembly/host_unmapped/*_scaffolds.fasta $outdir/05_Assembly/host_unmapped
rsync -r --progress -v $prefix/05_Assembly/host_unmapped/*_assembly_graph.fastg $outdir/05_Assembly/host_unmapped
rsync -r --progress -v $prefix/05_Assembly/host_unmapped/*_spades.log $outdir/05_Assembly/host_unmapped
# The headers will have an extra gnl|Prokka| prefixed to them but they can be parsed out with sed and related to contigs by looking into contig fates
rsync -r --progress -v $prefix/database/MAGs_database_bed_files $outdir/database
rsync -r --progress -v $prefix/database/MAGs_database_ffn_files $outdir/database
rsync -r --progress -v $prefix/database/MAGs_database_gff_files $outdir/database
rsync -r --progress -v $prefix/database/MAGs_database_faa_files $outdir/database
rsync -r --progress -v $prefix/database/MAGs_database_genomes $outdir/database
rsync -r --progress -v $prefix/database/MAGs_database_Orthofinder $outdir/database
