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
# rule make_MAG_reduced_db:
#     input:
#         mag_database = "database/MAGs_database",
#         ref_info = rules.select_cluster_ref_genomes_mag_database.output.ref_info,
#     output:
#         mag_database_reduced = "database/MAGs_database_reduced",
#         index = multiext("database/MAGs_database_reduced", ".amb", ".ann", ".bwt", ".pac", ".sa")
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = convertToMb("8G")
#     threads: 8
#     log: "logs/make_MAG_reduced_db.log"
#     benchmark: "logs/make_MAG_reduced_db.benchmark"
#     conda: "envs/core-cov-env.yaml"
#     shell:
#         """
#         python3 scripts/subset_metagenome_db.py --input {input.mag_database} --ref_info {input.ref_info}
#         bwa index {output.mag_database_reduced}
#         """
#
# rule subset_orthofiles_for_reduced_database:
#     input:
#         ref_info = rules.select_cluster_ref_genomes_mag_database.output.ref_info,
#         ortho_file = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho.txt",
#     output:
#         ortho_file = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho_reduced.txt",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = convertToMb("8G")
#     threads: 8
#     log: "logs/{group}_subset_orthofiles_for_reduced_database.log"
#     benchmark: "logs/{group}_subset_orthofiles_for_reduced_database.benchmark"
#     conda: "envs/core-cov-env.yaml"
#     shell:
#         """
#         python3 scripts/subset_orthofiles.py --input {input.orthofile} --ref_info {input.ref_info} --output {output.orthofile} --group {wildcards.group}
#         """
#
# rule map_to_MAGs_reduced:
#     input:
#         reads1 = rules.trim.output.reads1,
#         reads2 = rules.trim.output.reads2,
#         index = multiext("database/MAGs_database_reduced", ".amb", ".ann", ".bwt", ".pac", ".sa"),
#         genome_db = "database/MAGs_database_reduced"
#     output:
#         bam_mapped_reduced_temp = temp("09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}_mapped_temp.bam"),
#         bam_mapped_reduced = temp("09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}_mapped.bam"),
#         bam_mapped_reduced_index = temp("09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}_mapped.bam.bai"),
#         bam_reduced = temp("09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}.bam"),
#         sam_reduced = temp("09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}.sam"),
#         sam_mapped_reduced = temp("09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}_mapped.sam"),
#         mag_mapping_flagstat_reduced = "09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}_flagstat.tsv",
#     params:
#         perl_extract_mapped = "scripts/filter_sam_aln_length.pl",
#         mailto="aiswarya.prasad@unil.ch",
#         mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
#         jobname="{sample}_microbiomedb_mapping",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-4:10:00"),
#     resources:
#         mem_mb = convertToMb("20G")
#     threads: 8
#     log: "logs/{sample}_map_to_mags_reduced.log"
#     benchmark: "logs/{sample}_map_to_mags_reduced.benchmark"
#     conda: "envs/mapping-env.yaml"
#     shell:
#         """
#         bwa mem -t {threads} {input.genome_db} {input.reads1} {input.reads2} > {output.sam_reduced}
#         samtools view -bh {output.sam_reduced} | samtools sort - > {output.bam_reduced}
#         samtools flagstat -O tsv {output.bam_reduced} > {output.mag_mapping_flagstat_reduced}
#         samtools view -h -F4 -@ {threads} {output.sam_reduced} | samtools view - -F 0x800 -h | grep -E "NM:i:[0-4][[:blank:]]|^\@" | perl {params.perl_extract_mapped} - > {output.sam_mapped_reduced}
#         samtools view -bh {output.sam_mapped_reduced} | samtools sort - > {output.bam_mapped_reduced_temp}
#         picard -Xmx8g AddOrReplaceReadGroups I={output.bam_mapped_reduced_temp} O={output.bam_mapped_reduced} RGID={wildcards.sample} RGLB=lib1 RGPL=illumina RGPU=none RGSM={wildcards.sample}
#         samtools index {output.bam_mapped_reduced}
#         """
#
# rule de_duplicate:
#     input:
#         bam_red = rules.map_to_MAGs_reduced.output.bam_mapped_reduced,
#     output:
#         bam_red = temp("09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}_mapped_red_dedup.bam"),
#         bam_red_index = temp("09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}_mapped_red_dedup.bam.bai"),
#         txt = "09_MagDatabaseProfiling/SNVProfiling/Mapping/{sample}_dedup_metrics.txt"
#     params:
#         perl_extract_mapped = "scripts/filter_sam_aln_length.pl",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-6:10:00"),
#     resources:
#         mem_mb = 8000
#     threads:
#         20
#     conda: "envs/mapping-env.yaml"
#     log: "logs/{sample}_de_duplicate.log"
#     benchmark: "logs/{sample}_de_duplicate.benchmark"
#     shell:
#         """
#         picard -Xmx8g MarkDuplicates INPUT={input.bam_red} OUTPUT={output.bam_red} METRICS_FILE={output.txt};
#         samtools index {output.bam_red}
#         """
#
# rule core_cov_red:
#     input:
#         bed_files = lambda wildcards: expand("database/bed_files/{genome}.bed", genome=get_g_dict_for_groups_from_data_reduced_db(checkpoints.make_phylo_table.get().output.out_mags_filt, rules.select_cluster_ref_genomes_mag_database.output.ref_info)[wildcards.group]),
#         ortho_file = "database/MAGs_database_Orthofinder/{group}/OrthoFinder/{group}_single_ortho_reduced.txt",
#         bam_file = rules.de_duplicate.output.bam_red,
#         ref_info = rules.select_cluster_ref_genomes_mag_database.output.ref_info, # does not have to be subset because script only reads lines corresponsing to rep genomes
#     output:
#         core_cov_txt = "09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/{sample}/{group}_corecov.txt"
#     params:
#         bedfiles_dir = "database/bed_files/",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
#     log: "logs/{sample}_{group}_core_cov_red.log"
#     benchmark: "logs/{sample}_{group}_core_cov_red.benchmark"
#     conda: "envs/core-cov-env.yaml"
#     threads: 5
#     shell:
#         """
#         python3 scripts/core_cov.py --info {input.ref_info} --bamfile {input.bam_file} --ortho {input.ortho_file} --beddir {params.bedfiles_dir} --outfile {output.core_cov_txt} --group {wildcards.group} --sample {wildcards.sample}
#         """
#
# rule merge_core_cov_red:
#     input:
#         core_cov_per_sample = expand("09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/{sample}/{{group}}_corecov.txt", sample=SAMPLES)
#     output:
#         core_cov_txt = "09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/Merged/{group}_corecov.txt"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
#     log: "logs/{group}_merge_core_cov_red.log"
#     benchmark: "logs/{group}_merge_core_cov_red.benchmark"
#     run:
#         header_done = False
#         with open(output.core_cov_merged, "w") as out_fh:
#             for file in input.core_cov_per_sample:
#                 with open(file, "r") as in_fh:
#                     header = in_fh.readline()
#                     if not header_done:
#                         out_fh.write(header)
#                         header_done = True
#                     for line in in_fh:
#                         out_fh.write(line)
#
# rule core_cov_plots_red:
#     input:
#         txt = "09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/Merged/{group}_corecov.txt",
#     output:
#         txt = "09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/Merged/{group}_coord.txt"
#         pdf = "09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/Merged/{group}_coord.pdf"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-1:10:00"),
#     resources:
#         mem_mb = 8000
#     conda: "envs/core-cov-env.yaml"
#     log: "logs/{group}_core_cov_plots_red.log"
#     benchmark: "logs/{group}_core_cov_plots_red.benchmark"
#     threads: 1
#     shell:
#         """
#         ./scripts/core_cov.R {input.txt};
#         """

# rule corecov_split_by_magOTU:
#     input:
#         corecovs = lambda wildcards: ["09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/Merged/"+group+"_corecov.txt" for group in [x for x in get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt).keys() if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_mags_filt) > 2]],
#         coords = lambda wildcards: ["09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/Merged/"+group+"_coord.txt" for group in [x for x in get_g_dict_for_groups_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt).keys() if num_genomes_in_group(x, checkpoints.make_phylo_table.get().output.out_mags_filt) > 2]],
#     output:
#         magOTU_coords = expand("09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/Merged/{magOTU}_split_coord.txt", magOTU=get_g_magOTU_dict_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt, rules.select_cluster_ref_genomes_mag_database.output.ref_info).keys()),
#         magOTU_corecovs = expand("09_MagDatabaseProfiling/SNVProfiling/CoverageEstimation/Merged/{magOTU}_split_corecov.txt", magOTU=get_g_magOTU_dict_from_data(checkpoints.make_phylo_table.get().output.out_mags_filt, rules.select_cluster_ref_genomes_mag_database.output.ref_info).keys()),
#     log: "logs/corecov_split_by_magOTU.log"
#     threads: 1
#     run:
# rule generate_freebayes_regions:
#     input:
#         ref_idx = "database/MAGs_database_reduced",
#         index = "database/MAGs_database_reduced.fai",
#         bams = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam", sample=config["SAMPLES"]),
#         bam_files_index = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam.bai", sample=config["SAMPLES"]),
#         script = "scripts/fasta_generate_regions.py",
#     output:
#         regions = expand("database/freebayes_regions/"+DBs["microbiome_red"]+".{chrom}.region.{i}.bed", chrom=CHROMS, i = CHUNKS)
#     log:
#         "logs/generate_freebayes_regions.log"
#     params:
#         chunks = CHUNKS,
#         nchunks = N,
#         chroms = CHROMS, # names of fasta header in database
#         outdir = "database/freebayes_regions/"+DBs["microbiome_red"],
#     shell:
#         """
#         python3 {input.script} {input.index} {params.nchunks} --chunks --bed {params.outdir} --chromosome {params.chroms}
#         """
#
# rule freebayes_profiling:
#     input:
#         bam_files = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam", sample=config["SAMPLES"]),
#         bam_files_index = expand("07_SNVProfiling/Mapping_red/{sample}_mapped_red_dedup.bam.bai", sample=config["SAMPLES"]),
#         db_red_ind = "database/"+DBs["microbiome_red"]+".fai",
#         db_red = "database/"+DBs["microbiome_red"],
#         regions = "database/freebayes_regions/"+DBs["microbiome_red"]+".{chrom}.region.{i}.bed"
#     output:
#         vcf = temp("07_SNVProfiling/VCFs/{chrom}/variants.{i}.vcf")
#     params:
#         help_script = os.path.join(os.getcwd(), "scripts/fasta_generate_regions.py"), # if not in path
#         freebayes = "freebayes", # path
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("1-00:00:00"),
#     resources:
#         mem_mb = 10000
#     threads: 15
#     conda: "envs/snv-env.yaml"
#     log: "logs/freebayes_profiling_{chrom}_{i}.log"
#     shell:
#         """
#         freebayes --pooled-continuous --min-alternate-count 5 --min-alternate-fraction 0.1 --min-coverage 10 -f {input.db_red} -t {input.regions} {input.bam_files} > {output.vcf}
#         """
#
# rule concat_vcfs:
#     input:
#         calls = expand("07_SNVProfiling/VCFs/{chrom}/variants.{i}.vcf", chrom=CHROMS, i=CHUNKS)
#     output:
#         vcf = "07_SNVProfiling/freebayes_raw.vcf"
#     log:
#         "logs/concat_vcfs.log"
#     conda:
#         "envs/snv-env.yaml"
#     threads: 15
#     shell:
#         "bcftools concat {input.calls} | vcfuniq > {output.vcf}"
#
# rule vcf_filter_indels:
#     input:
#         vcf_raw = rules.concat_vcfs.output.vcf
#     output:
#         vcf_filt = "07_SNVProfiling/All_CandSNVs.vcf"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 8000
#     conda: "envs/snv-env.yaml"
#     log: "logs/vcf_filter_indels.log"
#     shell:
#         """
#         bcftools view --exclude-types indels {input.vcf_raw} > {output.vcf_filt}
#         """
#
# rule filt_bedfiles:
#     input:
#         metafile = "database/"+DBs["microbiome"]+"_metafile.txt",
#         sdp_coords = expand("07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_coord.txt", sdp=get_g_sdp_dict(DBs["microbiome"]).keys()),
#         sdp_corecovs = expand("07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_corecov.txt", sdp=get_g_sdp_dict(DBs["microbiome"]).keys()),
#         bedfile = "database/bed_files/{genome}.bed",
#         ortho_dir = "database/OrthoFiles_"+DBs["microbiome_red"]
#     output:
#         core_red_bed = "07_SNVProfiling/bed_files/{genome}_core_filt.bed"
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 2000
#     log: "logs/{genome}_filt_bedfiles.log"
#     script:
#         "scripts/filt_core_bed.py"
#
# rule calc_core_lengths:
#     input:
#         core_red_beds = expand("07_SNVProfiling/bed_files/{genome}_core_filt.bed", genome=get_g_list(DBs["microbiome_red"])),
#         db = rules.subset_ortho_and_db.output.db, # so that it is repeated if it changes - not used by script
#     output:
#         lengths = "07_SNVProfiling/table_core_length.txt"
#     params:
#         genomes_sdp_dict=get_g_sdp_dict(DBs["microbiome"]),
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     resources:
#         mem_mb = 2000
#     log: "logs/calc_core_length.log"
#     script:
#         "scripts/calc_filt_core_length.py"
#
# rule cat_bedfiles:
#     input:
#         bedfiles = expand("07_SNVProfiling/bed_files/{genome}_core_filt.bed", genome=get_g_list(DBs["microbiome_red"]))
#     output:
#         bed_all = "07_SNVProfiling/bed_files/all_filt.bed"
#     log: "logs/cat_bedfiles.log"
#     threads: 1
#     shell:
#         """
#         cat {input.bedfiles} > {output.bed_all}
#         """
#
# rule vcf_filter_intersect:
#     input:
#         vcf_filt = rules.vcf_filter_indels.output.vcf_filt,
#         bed_all = rules.cat_bedfiles.output.bed_all
#     output:
#         vcf_intersect = "07_SNVProfiling/Intersect_CandSNVs.vcf",
#     params:
#         vcfintersect = "vcfintersect",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     threads:
#         20
#     conda: "envs/snv-env.yaml"
#     log: "logs/vcf_filter_intersect.log"
#     shell:
#         """
#         {params.vcfintersect} --bed {input.bed_all} {input.vcf_filt} | vcfbreakmulti > {output.vcf_intersect}
#         """
#
# rule vcf_split_by_sdp:
#     input:
#         vcf = rules.vcf_filter_intersect.output.vcf_intersect,
#         metafile = "database/"+DBs["microbiome"]+"_metafile.txt",
#     output:
#         vcfs = expand("07_SNVProfiling/{sdp}_all_samples.vcf", sdp=get_g_sdp_dict(DBs["microbiome"]).keys())
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-0:30:00"),
#     threads:
#         4
#     log: "logs/vcf_split_by_sdp.log"
#     script:
#         "scripts/vcf_split_by_sdp.py"
#
# rule vcf_filter_samples:
#     input:
#         vcf = "07_SNVProfiling/{sdp}_all_samples.vcf",
#         coord = "07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_coord.txt",
#     output:
#         vcf = "07_SNVProfiling/{sdp}.vcf",
#     params:
#         script = "scripts/filt_vcf_samples.py",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     threads:
#         20
#     conda: "envs/snv-env.yaml"
#     log: "logs/{sdp}_vcf_filter_samples.log"
#     shell:
#         """
#         SAMPLES=$( python3 {params.script} {input.vcf} {input.coord} )
#         echo $SAMPLES
#         if [ -z "$SAMPLES" ]; then
#             echo "No samples in list!"
#             touch {output.vcf}
#         else
#             vcfkeepsamples {input.vcf} $SAMPLES > {output.vcf}
#         fi
#         """
#
# rule vcf_filter_missing:
#     input:
#         vcf = "07_SNVProfiling/{sdp}.vcf"
#     output:
#         freq = "07_SNVProfiling/{sdp}_filt.freq",
#     params:
#         summary = "07_SNVProfiling/summary_filtering.txt",
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-2:10:00"),
#     threads:
#         4
#     log: "logs/{sdp}_vcf_filter_missing.log"
#     script:
#         "scripts/filter_snvs.py"
#
# rule summarise_snps:
#     input:
#         metadata = "Metadata_"+ProjectIdentifier+".csv",
#         lengths = "07_SNVProfiling/table_core_length.txt",
#         freq = rules.vcf_filter_missing.output.freq,
#         coord = "07_SNVProfiling/CoreCov_"+ProjectIdentifier+"/{sdp}_split_coord.txt",
#     output:
#         txt_sample = "07_SNVProfiling/{sdp}_sample_var_host.txt",
#         txt_host = "07_SNVProfiling/{sdp}_tot_var.txt",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-0:20:00"),
#     threads:
#         4
#     conda: "envs/mapping-env.yaml" # for perl
#     log: "logs/{sdp}_summarise_snps.log"
#     script:
#         "scripts/summarise_snps.py"
#
# rule summarise_shared_fractions:
#     input:
#         freq = rules.vcf_filter_missing.output.freq
#     output:
#         shared_frac = "07_SNVProfiling/{sdp}_filt_shared_fraction.txt",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-3:10:00"),
#     threads:
#         4
#     log: "logs/{sdp}_summarise_shared_fractions.log"
#     conda: "envs/mapping-env.yaml" # for perl
#     shell:
#         """
#         perl scripts/calc_jaccard.pl {input.freq}
#         """
#
# rule make_cumulative_curves:
#     input:
#         metadata = "Metadata_"+ProjectIdentifier+".csv",
#         freq = rules.vcf_filter_missing.output.freq,
#         lengths = rules.calc_core_lengths.output.lengths,
#     output:
#         cum_curve = "07_SNVProfiling/{sdp}_cum_curve.txt",
#         cum_curve_by_group = "07_SNVProfiling/{sdp}_cum_curve_by_group.txt",
#         cum_curve_by_colony = "07_SNVProfiling/{sdp}_cum_curve_by_colony.txt",
#     params:
#         nb_curves = 10,
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-0:20:00"),
#     threads:
#         4
#     conda: "envs/mapping-env.yaml" # for perl
#     log: "logs/{sdp}_make_cumulative_curves.log"
#     script:
#         "scripts/cum_curve_SNVs_host.py"
#
# rule make_dist_matrix:
#     input:
#         file_frac = rules.summarise_shared_fractions.output.shared_frac
#     output:
#         dist_matrix = "07_SNVProfiling/{sdp}_dist_matrix.txt",
#     params:
#         mailto="aiswarya.prasad@unil.ch",
#         account="pengel_spirit",
#         runtime_s=convertToSec("0-0:20:00"),
#     threads:
#         4
#     conda: "envs/rmd-env.yaml"
#     log: "logs/{sdp}_make_dist_matrix.log"
#     script:
#         "scripts/distance_matrix.R"

checkpoint select_cluster_ref_genomes_mag_database:
    input:
        genomes_info = lambda wildcards: checkpoints.make_phylo_table.get().output.out_mags_filt,
        genome_scores = lambda wildcards: rules.drep.output.drep_S
    output:
        ref_info = "database/MAGs_database_ref_info.txt"
    params:
        mags_list = lambda wildcards: get_MAGs_list_dict(checkpoints.make_phylo_table.get().output.out_mags_filt, full_list = True),
        mailto="aiswarya.prasad@unil.ch",
        mailtype="BEGIN,END,FAIL,TIME_LIMIT_80",
        account="pengel_spirit",
        runtime_s=convertToSec("0-2:10:00"),
    benchmark: "logs/select_cluster_ref_genomes_mag_database.benchmark"
    log: "logs/select_cluster_ref_genomes_mag_database.log"
    run:
        MAGs_selected = set(params.mags_list)
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
                if genome_id not in MAGs_selected:
                    continue
                cluster = line.split("\t")[11]
                # genus = line.split("\t")[13]
                # genus = genus+cluster if genus == "g__" else genus
                if cluster in group_max_score.keys():
                    if int(genome_scores[genome_id]) > group_max_score[cluster]:
                        group_max_score[cluster] = genome_scores[genome_id]
                        rep_genome_dict[cluster] = genome_id
                else:
                    group_max_score[cluster] = genome_scores[genome_id]
                    rep_genome_dict[cluster] = genome_id
        with open(input.genomes_info, "r") as info_fh:
            with open(output.ref_info, "w") as out_fh:
                out_fh.write("genome_id\tcluster\trep_genome\tgroup\n")
                for line in info_fh:
                    if line.startswith("ID"):
                        continue
                    genome_id = line.split("\t")[0]
                    if "MAG_" not in genome_id:
                        continue
                    cluster = line.split("\t")[11]
                    genus = line.split("\t")[13]
                    if genus == "NA":
                        print(f"genus for {genome_id} is {genus}. This group may not make sense.")
                        continue
                    species_id = cluster
                    group_name = genus+cluster if genus == "g__" else genus
                    rep_genome = 1 if genome_id in rep_genome_dict.values() else 0
                    out_fh.write(f"{genome_id}\t{species_id}\t{rep_genome}\t{group_name}\n")
