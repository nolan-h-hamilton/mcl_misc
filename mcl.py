def cluster_features(mcl_output_file, species, normalize=False):

    """find pertinent features of clusters returned from mcl in one pass

    Args:
        mcl_output_file: txt file containing mcl clusters
        species: list of all species with sequences present in clusters
        normalize: True if species_frequency should be normalized
    Returns:
        A dictionary with keys:
            'species_frequency': a dict of species and with values denoting
                how often the respective species were present in clusters
                in `mcl_output_file`
            'cluster_sizes': a list containing sizes of all clusters ob-
                served
            'size_distr_dict': dict containing keys for mean, min, and
                 max of observed cluster sizes
            'singletons': a list of single-sequence clusters observed
            'num_clusters': total number of clusters observed
            'clusters': list of lists containing observed clusters
    Notes:
        This function should provide all data needed from the mcl file
            in a single traversal. If more features are of interest, this
            function should be modified to compute them without additional
            traversals (i.e., computed in the main loop)
    """
    clusters = []
    species_frequency = dict.fromkeys(species, 0)
    cluster_sizes = []
    singletons = []
    num_clusters = 0
    size_distr_dict = {'mean_clstr_size': 0.0, 'min_clstr_size': 0.0, 'max_clstr_size': 0.0}
    
    with open (mcl_output_file, 'r') as f:
        for line in f:
            # each 'lcl|' denotes cluster
            clusters = line.split('lcl|')
            for cluster in clusters:

                # get tab-separated seqs in cluster
                # and record cluster
                try:
                    cluster = cluster.split('\t')[0:-1]
                except AttributeError:
                    #blank line?
                    continue
                clusters.append(cluster)

                # this counter may not be needed, but may avoid potentially
                # expensive calls to len(clusters) in future. Keep for now.
                num_clusters += 1

                # compute cluster size features
                size = len(cluster)
                size_distr_dict['mean_clstr_size'] += size
                if size_distr_dict['min_clstr_size'] > size:
                    size_distr_dict['min_clstr_size'] = size
                if size_distr_dict['max_clstr_size'] < size:
                    size_distr_dict['max_clstr_size'] = size
                cluster_sizes.append(size)

                # record singleton clusters
                if size == 1:
                    singletons.append(cluster[0])

                # determine species of seqs in cluster
                # and increment freq for observed species
                present_specs = [x[0:7] for x in cluster]
                for spec in species:
                    if spec in present_specs:
                        species_frequency[spec] += 1

    # True,  divide frequency by num_clusters so that the value in dict
    # represents what fraction of clusters the species was present in
    if normalize:
        for key in species_frequency:
            species_frequency[key] /= num_clusters

    return {'species_frequency': species_frequency,
            'cluster_sizes': cluster_sizes,
            'size_distr_dict': size_distr_dict,
            'singletons': singletons,
            'num_clusters': num_clusters,
            'clusters': clusters}
