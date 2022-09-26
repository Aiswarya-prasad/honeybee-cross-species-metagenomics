#!/bin/bash

args=("$@")

scripts_dir=${args[0]}
genome_db_meta=${args[1]}
outfile=${args[2]}

for j in $( ls *.faa ); do
    OG=${j%."faa"}
    echo "Processing: $OG"
    aln_faa=$OG"_aln.fasta"
    #Aligning amino-acid sequences
    mafft --auto --quiet $j > $aln_faa
    #Back-translating alignment (codon-aligned nucleotide alignment)
    ffn_file=$OG".ffn"
    python3 "$scripts_dir/aln_aa_to_dna.py" "$aln_faa" "$ffn_file"
    #Trimming alignment for gaps
    aln_nuc=$OG"_aln_nuc.fasta"
    python3 "$scripts_dir/trim_aln.py"  "$aln_nuc"
    #Simplifying headers in alignment
    trim_file=$OG"_aln_trim.fasta"
    sed -i 's/_.*$//g' $trim_file
    #Calculating inter-SDP alignment stats
    python3 "$scripts_dir/calc_perc_id_orthologs.py" "$genome_db_meta" "$trim_file" "$outfile"
done
