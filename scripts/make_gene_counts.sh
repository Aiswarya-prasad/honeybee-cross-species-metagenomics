awk -F'\t' '

function counts_from_depths(depths_list, gene) {
    split(depths_list, depths_array)
    total_depth = 0
    for (i in depths_array) {
        total_depth += depths_array[i]
    }
    if (total_depth == 0) {
        return 0
    }
    zero_count = 0
    for (i in depths_array) {
        if (depths_array[i] == 0) {
            zero_count++
        }
    }
    if (zero_count > gene_lengths_dict[gene] * 0.5) {
        return 0
    }
    counts = median(depths_array)
    adjusted_count = counts / gene_lengths_dict[gene] * 1000
    return adjusted_count
}

BEGIN {
    OFS="\t"
}
NR == FNR {
    # Read gene_lengths_dict from file
    gene_lengths_dict[$1] = $2
    next
}
{
    gene = $1
    position = $2
    depth = $3
    gene_depths[gene] = gene_depths[gene] " " depth
}
END {
    for (gene in gene_depths) {
        print gene, counts_from_depths(gene_depths[gene], gene)
    }
}
' results/08_gene_content/03_gene_counts/gene_lengths.txt results/08_gene_content/01_profiling/F7-5_mapped.depth > results/08_gene_content/03_gene_counts/all_genes_matrix.txt