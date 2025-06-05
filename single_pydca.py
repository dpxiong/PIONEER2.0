import os, sys

def calculate_single_pydca(input_file, plmdca_output_dir, mfdca_output_dir, num_threads=6):
    os.system('plmdca compute_fn protein %s --max_iterations 500 --num_threads %s --output_dir %s --apc' % (input_file, str(num_threads), plmdca_output_dir))
    os.system('mfdca compute_fn protein %s --output_dir %s --apc' % (input_file, mfdca_output_dir))


if __name__ == '__main__':
    fasta_dict = {}
    with open('../data/all_preppi_uniprot_seqs.txt', 'r') as infile:
        for line in infile:
            line_list = line.strip().split()
            fasta_dict[line_list[0]] = line_list[1]
            
    with open('../data/for_sequence_extraction/updated_uniprots.txt', 'r') as infile:
        for line in infile:
            line_list = line.strip().split()
            assert line_list[0] in fasta_dict
            assert fasta_dict[line_list[0]] == line_list[1]
            
    print(len(fasta_dict))

    start = int(sys.argv[1])
    end = int(sys.argv[2])
    print(start, end)
    plmdca_output_dir = '../data/single_raw_plmdca'
    mfdca_output_dir = '../data/single_raw_mfdca'
    for identifier in sorted(list(fasta_dict.keys()))[start:end]:
        input_file = '../data/msa_uniref90/'+identifier+'.aligned_msa'
        if os.path.exists(input_file):
            print('>', identifier)
            calculate_single_pydca(input_file, plmdca_output_dir, mfdca_output_dir, 8)
        else:
            print('Not exist:', identifier)

    print('done!')
