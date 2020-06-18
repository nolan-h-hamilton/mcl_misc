def species_pwe_dict(abc_file, double_pairs=False, normalize=False):

    """form pairwise E-val dict of species in abc file

    Args:
        abc_file: tab-separated txt file of format SPECIES1 SPECIES2 E
        double_pairs: if True, keys are added for SPECIES1.SPECIES2 and
            SPECIES2.SPECIES1
        normalize: if True, the 0th index of each value in the list is
            divided by the 1st index to compute the average E value
    Returns:
        A dictionary of format {SPECIES1.SPECIES2: [SUM_E, SPEC_PAIR_CT]
            Note that the default beahvior is to create only one key for
            for each *unordered* pair of species, so  key 'AAMBTRI.ANUPAD
            might exist, but 'ANUPADV.AAMBTRI' will not.
    Notes:
        Can be run on abc files using B (bit-score) as well
    """

    pwe_dict = {}
    with open(abc_file, 'r') as f:
        for line in f:
            # process line text
            line = line.replace('lcl|', '').replace('\n','')
            split_line = line.split('\t')
            spec1 = get_tag(split_line[0])
            spec2 = get_tag(split_line[1])
            E = float(split_line[2])
            # make appropriate additions to pwe_dict
            # let's make one entry for each *unordered* pair
            if spec1 + '.' + spec2 not in pwe_dict:
                if spec2 + '.' + spec1 not in pwe_dict:
                    pwe_dict.update({spec1 + '.' + spec2: [E, 1]})
                else:
                    pwe_dict[spec2 + '.' + spec1][0] += E
                    pwe_dict[spec2 + '.' + spec1][1] += 1
            else:
                pwe_dict[spec1 + '.' + spec2][0] += E
                pwe_dict[spec1 + '.' + spec2][1] += 1

    if double_pairs:
        # .keys() returns an iterator, so cast to list or tuple
        for key in tuple(pwe_dict.keys()):
            split_key = key.split('.')
            pwe_dict.update({split_key[1] + '.' + split_key[0]: pwe_dict
                             [key]})
    if normalize:
        for key in pwe_dict:
            pwe_dict[key][0] /= pwe_dict[key][1]

    return pwe_dict
