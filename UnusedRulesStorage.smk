# targets:
    # midas_species_prevalence = "11_MIDAS/MidasProfiling/species_prevalence.txt",
    # midas_genes_summary = "11_MIDAS/midas_genes_summary.done",
    # snps_summary = "11_MIDAS/snps_summary.done",
    # midas_genus_species_prevalence = "11_MIDAS_genus_level/MidasProfiling/species_prevalence.txt",
    # midas_genus_genes_summary = "11_MIDAS_genus_level/midas_genes_summary.done",
    # snps_summary_genus = "11_MIDAS_genus_level/snps_summary.done",
    # snp_intra = lambda wildcards: expand("11_MIDAS/snp_diversity_intra/snp_diversity_{cluster}.info", cluster = get_cluster_dict(checkpoints.make_phylo_table.get().output.out_all).keys()),
    # snp_inter = lambda wildcards: expand("11_MIDAS/snp_diversity_inter/snp_diversity_{cluster}.info", cluster = get_cluster_dict(checkpoints.make_phylo_table.get().output.out_all).keys())

rule setup_midas_genus_parse_gff:
    input:
        gff_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.gff",
    output:
        outfile = "11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.genes",
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    benchmark: "logs/setup_midas_genus_parse_gff_{genome}.benchmark"
    log: "logs/setup_midas_genus_parse_gff_{genome}.log"
    shell:
        """
        python3 scripts/parse_gff_for_midas.py --outfile {output.outfile} --gff {input.gff_file}
        """

rule setup_midas_genus_for_custom_db_prepare_files:
    input:
        fna_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.fna",
        faa_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.faa",
        ffn_file = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.ffn",
    output:
        fna_file = "11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.fna",
        faa_file = "11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.faa",
        ffn_file = "11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.ffn"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    log: "logs/setup_midas_genus_for_custom_db_prepare_files_{genome}.log"
    shell:
        """
        mkdir -p 11_MIDAS_genus_level/MAGs_database/
        mkdir -p 11_MIDAS_genus_level/MAGs_database/{wildcards.genome}
        rsync -av {input.fna_file} {output.fna_file}
        rsync -av {input.faa_file} {output.faa_file}
        rsync -av {input.ffn_file} {output.ffn_file}
        """

rule setup_midas_genus_for_custom_db_write_mapfile:
    input:
        genomes_info = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
        genome_scores = lambda wildcards: rules.drep.output.drep_S,
    output:
        mapfile = "11_MIDAS_genus_level/midas_mapfile.tsv"
    params:
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    benchmark: "logs/setup_midas_genus_for_custom_db_write_mapfile.benchmark"
    log: "logs/setup_midas_genus_for_custom_db_write_mapfile.log"
    run:
        genome_scores = {}
        with open(input.genome_scores, "r") as win_fh:
            for line in win_fh:
                line = line.strip()
                if line.startswith("genome"):
                    continue
                genome_scores[line.split(",")[0].split(".fa")[0]] = float(line.split(",")[1])
        with open(input.genomes_info, "r") as info_fh:
            # make rep genomes set
            group_max_score = {}
            rep_genome_dict = {}
            for line in info_fh:
                if line.startswith("ID"):
                    continue
                genome_id = line.split("\t")[0]
                if "MAG_" not in genome_id:
                    continue
                cluster = line.split("\t")[11]
                genus = line.split("\t")[13]
                genus = genus+cluster if genus == "g__" else genus
                if genus in group_max_score.keys():
                    if int(genome_scores[genome_id]) > group_max_score[genus]:
                        group_max_score[genus] = genome_scores[genome_id]
                        rep_genome_dict[genus] = genome_id
                else:
                    group_max_score[genus] = genome_scores[genome_id]
                    rep_genome_dict[genus] = genome_id
        with open(input.genomes_info, "r") as info_fh:
            with open(output.mapfile, "w") as out_fh:
                out_fh.write("genome_id\tspecies_id\trep_genome\n")
                for line in info_fh:
                    if line.startswith("ID"):
                        continue
                    genome_id = line.split("\t")[0]
                    if "MAG_" not in genome_id:
                        continue
                    cluster = line.split("\t")[11]
                    genus = line.split("\t")[13]
                    # if genus == "NA":
                    #     print(f"genus for {genome_id} is {genus}.")
                    #     continue
                    species_id = genus+cluster if genus == "g__" else genus
                    rep_genome = 1 if genome_id in rep_genome_dict.values() else 0
                    out_fh.write(f"{genome_id}\t{species_id}\t{rep_genome}\n")

rule create_midas_genus_for_custom_db:
    input:
        fna_files = lambda wildcards: expand("11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.fna", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        faa_files = lambda wildcards: expand("11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.faa", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        ffn_files = lambda wildcards: expand("11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.ffn", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        genes_file = lambda wildcards: expand("11_MIDAS_genus_level/MAGs_database/{genome}/{genome}.genes", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
        mapfile = "11_MIDAS_genus_level/midas_mapfile.tsv"
    output:
        done = touch("11_MIDAS_genus_level/Midas_DB.done")
    params:
        outdir = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        # consider putting these paths in config
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        # midas_db_default = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/midas_db_v1.2/",
        indir = "11_MIDAS_genus_level/MAGs_database/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-23:10:00"),
    benchmark: "logs/create_midas_genus_for_custom_db.benchmark"
    log: "logs/create_midas_genus_for_custom_db.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        build_midas_db.py {params.indir} {input.mapfile} {params.outdir} --threads {threads} --resume
        if [ ! -f {params.midas_db}/marker_genes/phyeco.fa.bwt ];
        then
            hs-blastn index {params.midas_db}/marker_genes/phyeco.fa
        fi
        """

rule midas_genus_species_abundance:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        done = "11_MIDAS_genus_level/Midas_DB.done"
    output:
        output = "11_MIDAS_genus_level/MidasProfiling/{sample}/species/species_profile.txt",
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_species_abundance_{sample}.benchmark"
    log: "logs/midas_genus_species_abundance_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        run_midas.py species {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_genus_species_abundance_merge:
    input:
        per_sample_output = expand("11_MIDAS_genus_level/MidasProfiling/{sample}/species/species_profile.txt", sample=SAMPLES)
    output:
        relative_abundance = "11_MIDAS_genus_level/MidasProfiling/relative_abundance.txt",
        count_reads = "11_MIDAS_genus_level/MidasProfiling/count_reads.txt",
        coverage = "11_MIDAS_genus_level/MidasProfiling/coverage.txt",
        species_prevalence = "11_MIDAS_genus_level/MidasProfiling/species_prevalence.txt",
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_species_abundance_merged.benchmark"
    log: "logs/midas_genus_species_abundance_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py species {params.outdir} -i {params.outdir} -t dir
        """

rule midas_genus_gene_content:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        output = "11_MIDAS_genus_level/MidasProfiling/{sample}/species/species_profile.txt",
        done = "11_MIDAS_genus_level/Midas_DB.done",
    output:
        output = "11_MIDAS_genus_level/MidasProfiling/{sample}/genes/summary.txt",
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        species_cov = 10, # default is 3 use later if needed
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_gene_content_{sample}.benchmark"
    log: "logs/midas_genus_gene_content_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        run_midas.py genes {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_genus_gene_content_merge:
    input:
        per_sample_output = expand("11_MIDAS_genus_level/MidasProfiling/{sample}/genes/summary.txt", sample=SAMPLES)
    output:
        # genes_coverage = "11_MIDAS_genus_level/MidasProfiling/{cluster}/genes_coverage.txt",
        # genes_presence = "11_MIDAS_genus_level/MidasProfiling/{cluster}/genes_presence_absence.txt",
        # genes_summary = "11_MIDAS_genus_level/MidasProfiling/{cluster}/genes_summary.txt",
        # genes_copy_number = "11_MIDAS_genus_level/MidasProfiling/{cluster}/genes_copy_number.txt"
        outdone = touch("11_MIDAS_genus_level/midas_genes_summary.done")
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_gene_content_merged.benchmark"
    log: "logs/midas_genus_gene_content_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py genes {params.outdir} -i {params.outdir} -t dir
        for readme in 11_MIDAS_genus_level/MidasProfiling/*/readme.txt;
        do
            cp ${{readme%%readme.txt}}/readme.txt ${{readme%%readme.txt}}/genes_readme.txt
        done
        """

rule midas_genus_profile_snps:
    input:
        reads1 = rules.host_mapping_extract_host_filtered_reads.output.reads1,
        reads2 = rules.host_mapping_extract_host_filtered_reads.output.reads2,
        output = "11_MIDAS_genus_level/MidasProfiling/{sample}/species/species_profile.txt",
        done = "11_MIDAS_genus_level/Midas_DB.done"
    output:
        summary = "11_MIDAS_genus_level/MidasProfiling/{sample}/snps/summary.txt",
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/{sample}/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        species_cov = 1, # default is 3 use later if needed
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_profile_snps_{sample}.benchmark"
    log: "logs/midas_genus_profile_snps_{sample}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        run_midas.py snps {params.outdir} -1 {input.reads1} -2 {input.reads2} -t {threads} -d {params.midas_db} --remove_temp
        """

rule midas_genus_profile_snps_merge:
    input:
        per_sample_output = expand("11_MIDAS_genus_level/MidasProfiling/{sample}/snps/summary.txt", sample=SAMPLES)
    output:
        # snps_freq = "11_MIDAS_genus_level/MidasProfiling/{cluster}/snps_freq.txt",
        # snps_depth = "11_MIDAS_genus_level/MidasProfiling/{cluster}/snps_depth.txt",
        # snps_summary = "11_MIDAS_genus_level/MidasProfiling/{cluster}/snps_summary.txt",
        # snps_info = "11_MIDAS_genus_level/MidasProfiling/{cluster}/snps_info.txt"
        outdone = touch("11_MIDAS_genus_level/snps_summary.done")
    params:
        outdir = "11_MIDAS_genus_level/MidasProfiling/",
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_genus_profile_snps_merged.benchmark"
    log: "logs/midas_genus_profile_snps_merged.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        merge_midas.py snps {params.outdir} -i {params.outdir} -t dir
        for readme in 11_MIDAS_genus_level/MidasProfiling/*/readme.txt;
        do
            cp ${{readme%%readme.txt}}/readme.txt ${{readme%%readme.txt}}/snps_readme.txt
        done
        """

rule midas_genus_snp_diversity_intra:
    input:
        done = rules.midas_profile_snps_merge.output.outdone
    output:
        outfile = "11_MIDAS_genus_level/snp_diversity_intra/snp_diversity_{cluster}.info"
    params:
        indir = "11_MIDAS_genus_level/MidasProfiling/{cluster}/",
        maf = 0.01, #minimum minor allele frequency to be considered polymorphic
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_snp_diversity_intra_{cluster}.benchmark"
    log: "logs/midas_snp_diversity_intra_{cluster}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        snp_diversity.py {params.indir} --genomic_type genome-wide --sample_type per-sample --out {output.outfile} --snp_maf {params.maf}
        """

rule midas_genus_snp_diversity_inter:
    input:
        done = rules.midas_profile_snps_merge.output.outdone
    output:
        outfile = "11_MIDAS_genus_level/snp_diversity_inter/snp_diversity_{cluster}.info"
    params:
        indir = "11_MIDAS_genus_level/MidasProfiling/{cluster}/",
        maf = 0.01, #minimum minor allele frequency to be considered polymorphic
        midas_db = "11_MIDAS_genus_level/MAGs_database/Midas_DB",
        midas_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/",
        midas_scripts_path = "/work/FAC/FBM/DMF/pengel/spirit/aprasad/Software/MIDAS/scripts/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-5:10:00"),
    benchmark: "logs/midas_snp_diversity_inter_{cluster}.benchmark"
    log: "logs/midas_snp_diversity_inter_{cluster}.log"
    threads: 8
    resources:
        mem_mb = 8000
    conda: "envs/midas-env.yaml"
    shell:
        """
        export PYTHONPATH={params.midas_path}
        export PATH=$PATH:{params.midas_scripts_path}
        export MIDAS_DB={params.midas_db}
        snp_diversity.py {params.indir} --genomic_type genome-wide --sample_type pooled-samples --out {output.outfile} --snp_maf {params.maf}
        """


rule prepare_genomes_for_metapop:
    input:
        mag = "07_AnnotationAndPhylogenies/01_prokka/{genome}/{genome}.fna",
    output:
        mag = "database/MAGs_database_genomes/{genome}.fna"
    params:
        mailto="aiswarya.prasad@unil.ch",
        account="pengel_spirit",
        runtime_s=convertToSec("0-1:10:00"),
    resources:
        mem_mb = 8000
    log: "logs/{genome}_prepare_genomes_for_metapop.log"
    threads: 4
    shell:
        """
        cat {input.mag} > {output.mag}
        """

rule make_norm_file_metapop:
    input:
        flagstats = expand("09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_flagstat.tsv", sample=SAMPLES)
    output:
        norm = "10_metapop/norm_file.tsv"
    log: "logs/make_norm_file_metapop.log"
    run:
        with open(output.norm, "w") as out_fh:
            for flagstat in input.flagstats:
                sample = flagstat.split("/")[-1].split("_flagstat.tsv")[0]
                with open(flagstat, "r") as flagstat_fh:
                    for line in flagstat_fh:
                        line = line.strip()
                        line_split = line.split("\t")
                        if line_split[2] == "primary":
                            reads = line_split[0]
                    out_fh.write(f"{sample}\t{reads}\n")

rule prepare_bams_metapop:
    input:
        bam_file = "09_MagDatabaseProfiling/MAGsDatabaseMapping/{sample}_mapped.bam",
    output:
        # move away file and touch to trick snakemake
        bam_file = "10_metapop/Bamfiles/{sample}.bam"
    params:
        bam_files_dir = "10_metapop/Bamfiles/",
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("1-4:10:00"),
    resources:
        mem_mb = convertToMb("4G")
    threads: 4
    log: "logs/{sample}_prepare_bams_metapop.log"
    shell:
        """
        mkdir -p {params.bam_files_dir}
        cp -n {input.bam_file} {output.bam_file}
        """

rule run_metapop:
    input:
        bam_files = expand("10_metapop/Bamfiles/{sample}.bam", sample=SAMPLES),
        norm = "10_metapop/norm_file.tsv",
        reference_genomes = lambda wildcards: expand("database/MAGs_database_genomes/{genome}.fna", genome=get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True)),
    output:
        out = touch("10_metapop/metapop.done")
    params:
        prefix = os.getcwd(),
        outdir = "10_metapop/",
        reference_genomes_dir = "database/MAGs_database_genomes/",
        bam_files_dir = "10_metapop/Bamfiles/",
        aln_length = 50, # default is 30
        min_cov = 10, # contigs with breadth of coverage (#bases covered/contig length) less than this are removed from microdiversity.  default 20
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("1-4:10:00"),
    resources:
        mem_mb = convertToMb("50G")
    threads: 16
    log: "logs/run_metapop.log"
    benchmark: "logs/run_metapop.benchmark"
    conda: "envs/snv-env.yaml"
    shell:
        """
        #
        # some lines to run because metapop fails to do this for some reason
        mkdir -p 10_metapop
        mkdir -p 10_metapop/MetaPop
        mkdir -p 10_metapop/MetaPop/01.Genomes_and_Genes
        # ran on the from end
        # awk \'/^>/ {{printf(\"\\n%s\\n\",$0);next; }} {{ printf(\"%s\",$0);}}  END {{printf(\"\\n\");}}\' < database/MAGs_database | tail -n +2 > 10_metapop/MetaPop/01.Genomes_and_Genes/all_genomes.fasta
        prodigal -p meta -i 10_metapop/MetaPop/01.Genomes_and_Genes/all_genomes.fasta -d 10_metapop/MetaPop/01.Genomes_and_Genes/all_genomes_genes.fasta -o 10_metapop/MetaPop/01.Genomes_and_Genes/temp.txt
        #
        metapop --input_samples {params.prefix}/{params.bam_files_dir} --reference {params.prefix}/{params.reference_genomes_dir} --norm {params.prefix}/{input.norm} --threads {threads} --min_len {params.aln_length} --library $CONDA_PREFIX/lib/R/library --output {params.outdir}
        """
