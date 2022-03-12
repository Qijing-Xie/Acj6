"""
This file reads in a Acj6 binding site consensus file
and a txt file containing promoters of candidate genes,
calculate if promoters of those genes have Acj6 binding site(s),
and outputs the putative binding site in a .csv file
"""
import pandas as pd
from Bio.Seq import Seq
import sys

def get_gene_name(nameline):
    """
    input: a line contianing gene name and promoter site info
    (https://epd.epfl.ch/drosophila/drosophila_database.php?db=drosophila)
    output: name of gene
    """
    start = nameline.find(' ') + 1
    end = nameline.find(' ', start+1)

    name = nameline[start:end]
    return name


def get_gene_promoters(filename):
    """
    input: filename of file containing the promoter region
    output: dictionary of gene names as keys and the promoter region of genes as values
    """
    gene_names = []
    promoter_seqs = []
    file = open(filename,mode='r')
    i = 0
    for line in file:
        Lstripped=line.strip().upper()
        # check to see if the line represent gene name
        if Lstripped[0] == '>':
            i += 1  # change the indentation for gene sequence input
            gene_names.append(get_gene_name(line))
        else:
            # for same gene, stich sequence together
            if i > len(promoter_seqs):
                promoter_seqs.append(Lstripped)
            else:
                promoter_seqs[i-1] = promoter_seqs[i-1] + Lstripped
    file.close()
    dict_genename_promoters = dict(zip(gene_names, promoter_seqs))
    return dict_genename_promoters


def get_consensus_score(seq, acj6_consensus):
    """
    input: a 12bp sequence and acj6_binding site
    output: a score for Acj6 binding site
    """
    score = 0
    # calculate score for current 12 bp sequence
    for pos in acj6_consensus.columns:
        nucleotide = seq[acj6_consensus.columns[pos-1]-1] # get nucleotide at that position
        score += acj6_consensus.loc[nucleotide, pos]

    return score


def get_all_binding_sites(acj6_consensus, dict_genename_promoters):
    """
    input: consensus for Acj6, and dictionary containing genenames and promoters
    return: pandas dataframe containing all significant binding sites
    """
    df_binding_sites = pd.DataFrame(columns=['gene_name', 'Acj6_site', 'score'])
    # iterate through all genes in list
    binding_site_length = acj6_consensus.shape[1]
    threshold = 6 # based on Bai et al. 2009
    for gene in dict_genename_promoters.keys():
        promoter = dict_genename_promoters[gene] # get promoter sequence
        scan_len = len(promoter) - binding_site_length + 1  ## find max scanning length

        k = 0 ## initialize count for all possible binding sequences on promoter
        while k < scan_len: # loop through every possible binding sequence
            seq = promoter[k:k+binding_site_length]
            anti_seq = str(Seq(seq).reverse_complement())
            score = get_consensus_score(seq, acj6_consensus)
            score_anti = get_consensus_score(anti_seq, acj6_consensus)

            # add sequence information to df_binding_site if pass threshold
            if score > threshold:
                df_binding_sites = df_binding_sites.append({'gene_name': gene,
                                                             'Acj6_site': seq,
                                                             'score': score}
                                                            , ignore_index=True)
            elif score_anti > threshold:
                df_binding_sites = df_binding_sites.append({'gene_name': gene,
                                                             'Acj6_site': anti_seq,
                                                             'score': score_anti}
                                                            , ignore_index=True)
            k += 1
    return df_binding_sites


def main():
    args = sys.argv[1:]

    # get file containing promoter and pd dataframe containing the acj6 consnsus
    file_name = args[0]
    acj6_consensus = pd.read_csv(args[1], index_col=0)
    acj6_consensus.columns = acj6_consensus.columns.astype(int) # change column name to int

    #
    dict_genename_promoters = get_gene_promoters('input/'+file_name)
    df_sites = get_all_binding_sites(acj6_consensus, dict_genename_promoters)

    # save all binding sites
    output_filename = 'output/' + file_name[:file_name.find('.')] + '_out.csv'
    df_sites.to_csv(output_filename)

if __name__ == '__main__':
    main()
