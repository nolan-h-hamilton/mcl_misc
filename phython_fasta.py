"""Some functions for working with FASTA files"""

def get_tag(seq):
    """get species tag from a sequence name

    Args:
        seq: sequence name delimited by '.' and/or '_' with species tag in pos 0

    Raises:
        ValueError: if seq is None or does not contain '.' or '_'

    Returns:
        The species tag of sequence name `seq`

    Note:
        this function looks for the *first* delimiter and
        returns a string of characters preceding it. This
        method of getting species tags will prevent most
        bugs caused by non-uniform sequence name formatting
    """

    pos_first_delim = -1
    # consider adding 0-9 to `delimiters` so support TAG[0-9]
    # style formatting of sequence names
    delimiters = ['|', '_', '.']
    for i, char in enumerate(seq):
        if char in delimiters:
            pos_first_delim = i
            break
        
    if pos_first_delim == -1:
        raise ValueError('get_tag(): sequence names must be delimited by "." or "_"')

    return seq[0:pos_first_delim]


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

