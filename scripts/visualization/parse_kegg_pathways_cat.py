import os
import pandas as pd

main_sections = ['Metabolism',
                 'Genetic Information Processing',
                 'Environmental Information Processing',
                 'Cellular Processes',
                 'Organismal Systems',
                 'Human Diseases',
                 'Drug Development']

def clean_path_name(path_name_str):
    if path_name_str.startswith('M R '):
        path_name_str = path_name_str.split('M R ')[1]
    else:
        if path_name_str.startswith('M N '):
            path_name_str = path_name_str.split('M N ')[1]
        else:
            if path_name_str.startswith('M '):
                path_name_str = path_name_str.split('M ')[1]
            if path_name_str.startswith('N '):
                path_name_str = path_name_str.split('N ')[1]
            if path_name_str.startswith('R '):
                path_name_str = path_name_str.split('R ')[1]
    return path_name_str

# copied from https://www.genome.jp/kegg/pathway.html and modified by hand
with open('scripts/visualization/kegg_pathway_categories.tsv', 'w+') as f_out:
    f_out.write(f'section\tsubsection\tpath_name\n')
    path_name = ''
    subsection = ''
    section = ''
    with open('scripts/visualization/kegg_pathway_subcategories.txt', 'r') as f:
        for line in f:
            if line.startswith('\t\t'):
                path_name = ' '.join(line.strip().split()[1:])
                path_name = clean_path_name(path_name)
            else:
                if line.startswith('\t'):
                    subsection = line.strip()
                else:
                    section = line.strip()
            if subsection == '' or section == '' or path_name == '':
                continue
            else:
                f_out.write(f'{section}\t{subsection}\t{path_name}\n')
