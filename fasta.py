
def fasta_dict(fasta_file):
    fa_dict = {}
    with open(fasta_file, 'r') as fa_file:
        in_seq_data = False
        seq_name = None
        seq_data = None
        for line in fa_file:
            line = line.strip()
            if line[0] == '>':
                in_seq_data = True
                seq_name = line[1::]
                continue
            if in_seq_data:
                seq_data = line
                fa_dict.update({seq_name: seq_data})
                in_seq_data = False
    return fa_dict

