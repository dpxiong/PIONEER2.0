if __name__ == '__main__':
    import os, sys
    from msa import *
    
    fasta_dict = {}
    with open('../data/all_preppi_uniprot_seqs.txt', 'r') as infile:
        for line in infile:
            line_list = line.strip().split()
            fasta_dict[line_list[0]] = line_list[1]
    
    interactions = []
    with open('../data/all_preppi_interactions.txt', 'r') as infile:
        for line in infile:
            line_list = line.strip().split()
            interactions.append((line_list[0], line_list[1]))
            
    assert len(interactions) == len(set(interactions))
    
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    print(start, end)
    for interaction in interactions[start:end]:
        p1, p2 = interaction[0], interaction[1]
        print('>', p1, p2)
        cdhit_input_file1 = os.path.join('../data/msa', p1+'_rawmsa.fasta')
        cdhit_input_file2 = os.path.join('../data/msa', p2+'_rawmsa.fasta')
        cdhit_clstr_file1 = os.path.join('../data/msa', p1+'.cdhit.clstr')
        cdhit_clstr_file2 = os.path.join('../data/msa', p2+'.cdhit.clstr')
        if os.path.exists(cdhit_input_file1) and os.path.exists(cdhit_input_file2) and os.path.exists(cdhit_clstr_file1) and os.path.exists(cdhit_clstr_file2):
            joined_oneline_msa_file, joined_aligned_msa_file = calculate_joined_msa(p1, p2, fasta_dict[p1], fasta_dict[p2], cdhit_input_file1, cdhit_input_file2, cdhit_clstr_file1, cdhit_clstr_file2, '../data/joined_msa')
        
    print('done!')
