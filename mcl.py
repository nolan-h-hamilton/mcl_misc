
"""Compute some basic features of MCL results

    Description:
        cluster_features() is responsible for the majority of functionality.
        pertinent features/data from the mcl results located in the dump file
        generated via mclblastline, mcxdump, etc. are traversed and computed in
        a single-pass whenever possible and returned as a dictionary.
         
    Some Vocabulary:

        -bijective cluster: clusters which contain len(species)
             sequences and which have a sequence from each species
        
        -complete cluster: a cluster which contains at least one sequence
             for every provided species (args.species)

        -singleton cluster: a cluster containing only one sequence

        -excluded species: a species whose exclusion in a particular cluster
             prevents it from being a bijection. for example,
             if species = "ANUPADV, AAMBTRI, ALIRTUL" and cluster = "ANUPADV.1001,
             AAMBTRI.2833", then "ALIRTUL" is the "excluded" species. May want
             to come up with less ambiguous term for this...

        -tag: the tag name of a species used to format sequence names
            'AAMBTRI' for the sequence of Amborella trichopoda, 'AAMBTRI.1001'
        
    Notes:
        -mcl results are assumed to be in the format of an mcxdump dump file.
             this file is generated via "mclblastline" and is specified with
             prefix "dump". (https://micans.org/mcl/man/mclblastline.html)

        -species tags will be interpreted as the string of characters preceding
             the first delimiter in sequence names. For example, for sequence
             name "ANUPADV.1001", the species tag will be interpreted as "ANUPADV".
             supported delimiters are ".", "_", "|". Multiple delimiters can be
             used in sequence name formats to specify other aspects of seqs, but
             the first delimiter will mark the end of the species tag. That is,
                 "AAMB.00061_54" --> AAMB
                 "ANU|00061.54" --> ANU
                 "FLYGJAP_0061.64" --> FLYGJAP

        -for more info, see mcl_misc/data and/or the cluster_features()
             docstring


    Basic Usage Example:
        # in the following examples, mcl_sample.out is an mcl dumpfile
        # and tags_sample.txt is a text file containing chosen species
        # tags separated on by newlines. See mcl_misc/data/ for examples.
         
        # generate mcl dumpfile features with cluster_features()
        python mcl.py data/mcl_sample.out --species data/tags_sample.txt

        # add "--plot" to generate a log-scaled histogram of cluster sizes
        # and other features

        # add "--text" to print mcl result features to terminal

        # add "--print_bijections" to print clusters satisfying
        # the criteria of a bijection as defined above

        # add "--fasta_file [PATH_TO_FASTA] to generate an MFASTA
        # file from the results using "complete" clusters.

        # add "--normalize" to normalize several features by 
        # total number of clusters (e.g., species_freq_dict)


        # all of the parameters can be combined for complete output
"""


import sys
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import fasta


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
            'num_seqs': number of observed unique sequences
            'clusters': list of lists containing observed clusters
            'species_bijections': a list containing clusters with len(species) sequences with 
                 every species represented--clusters that exhibit a one-to-one mapping between
                 sequences and species
            'completes': a list containing clusters which contain sequences from all specified species
            'excluded_species_cts': a dictionary containing species keys and how often the species'
                 absence from a cluster was responsible for the cluster's disqualification as a "bijection"
             
    Notes:
        This function should provide all data needed from the mcl file
            in a single traversal. If more features are needed, this
            function should be modified to compute them without additional
            traversals (i.e., computed in the main loop) if possible.
    """
    
    clusters = []
    species_frequency = dict.fromkeys(species, 0)
    excluded_species_cts = dict.fromkeys(species, 0)
    cluster_sizes = []
    species_bijections = []
    singletons = []
    completes= []
    num_clusters = 0
    num_seqs = 0
    size_distr_dict = {'mean_clstr_size': 0.0, 'min_clstr_size': 100000.0,
                       'max_clstr_size': 0.0}
    
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
                    continue

                size = len(cluster)
                if size == 0:
                    continue
                
                # this count will assume a seq cannot be present in multiple clusters
                num_seqs += size
                
                # this counter may not be needed, but may avoid potentially
                # expensive calls to len(clusters) in future. Keep for now.
                num_clusters += 1
                clusters.append(cluster)

                size_distr_dict['mean_clstr_size'] += size
                if size_distr_dict['min_clstr_size'] > size:
                    size_distr_dict['min_clstr_size'] = size
                if size_distr_dict['max_clstr_size'] < size:
                    size_distr_dict['max_clstr_size'] = size
                cluster_sizes.append(size)

                # record singleton clusters
                if size == 1:
                    singletons.append(cluster[0])

                # compute species representation features
                # as of 06/17/2020, species tags can be any
                # length, but must not contain any delimiters
                present_specs = [get_tag(seq) for seq in cluster]
                set_species = set(species)
                set_pres = set(present_specs)
                spec_diff = set_species - set_pres

                # record clusters that contain all species in `completes`
                if len(spec_diff) == 0:
                    completes.append(cluster)
                    
                # if there are len(species) sequences present and every
                # species is represented, we define the cluster to be a
                # "bijection" (a 1-to-1 mapping from species to sequences)
                if len(spec_diff) == 0 and len(cluster) == len(species):
                    species_bijections.append(cluster)
                    
                # if there is only one species excluded that is preventing
                # the cluster from being a bijection, record the species
                elif len(spec_diff) == 1 and len(cluster) == len(species)-1:
                    excluded_species_cts[spec_diff.pop()] += 1

                #increment frequency for each present species
                for spec in set_pres:
                    species_frequency[spec] += 1


    size_distr_dict['mean_clstr_size'] /= num_clusters
    
    # True,  divide frequency by num_clusters so that the value in dict
    # represents what fraction of clusters the species was present in
    if normalize:
        for key in species_frequency:
            species_frequency[key] /= num_clusters

        off_by_one_ct = float(sum(excluded_species_cts.values()))
        for key in excluded_species_cts:
            excluded_species_cts[key] /= off_by_one_ct
            
        # if the user has opted to normalize, add a key for the number
        # of "off_by_one" clusters
        excluded_species_cts.update({'off_by_one_ct': off_by_one_ct})

    return {'species_frequency': species_frequency,
            'cluster_sizes': cluster_sizes,
            'size_distr_dict': size_distr_dict,
            'singletons': singletons,
            'completes': completes,
            'num_clusters': num_clusters,
            'num_seqs': num_seqs,
            'species_bijections': species_bijections,
            'excluded_species_cts': excluded_species_cts,
            'clusters': clusters}


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('mcl_file', help='file containing MCL results')
    parser.add_argument('--species', help='file containing species tags')
    parser.add_argument('--plot', help='display histogram of clusters features',
                        action='store_true')
    parser.add_argument('--text', help='print cluster features to terminal',
                        action='store_true')
    parser.add_argument('--normalize', help='normalize species frequencies',
                        action='store_true')
    parser.add_argument('--print_bijections', help='print a list of clusters that are one-to-one-mappings',
                        action='store_true')
    parser.add_argument('--print_completes', help='print clusters that contain all species',
                        action='store_true')
    parser.add_argument('--fasta_file', default=None,
                        help='set if user wishes to create mfasta file from mcl dump')
    
    args = parser.parse_args()

    species = [tag.strip() for tag in open(args.species, 'r').readlines()]
                        
    features_dict = cluster_features(args.mcl_file, species,
                                      normalize=args.normalize)

    clstr_size_arr = np.array(features_dict['cluster_sizes'])
    mean = round(features_dict['size_distr_dict']['mean_clstr_size'],2)
    max_ = round(features_dict['size_distr_dict']['max_clstr_size'],2)
    min_ = round(features_dict['size_distr_dict']['min_clstr_size'],2)
    singleton_prop = round(len(features_dict['singletons'])
                           / float(features_dict['num_clusters']),2)
    
    # these features can be computed/estimated in a single-pass
    # should consider adding iterative, incremental steps to 
    # to cluster_features() loop.
    std_dev = round(np.std(clstr_size_arr),2)
    med = round(np.median(clstr_size_arr),2)

    if args.plot:
        # SET Y-AXIS SCALE (linear, log, symlog, logit)
        plt.yscale('log')

        # SET X-AXIS RANGE
        x_range = 3*len(species) + 1
        plt.xlim(0, x_range)
        
        ylabel_string = ('FREQ\n\nmin: {}\nmedian: {}\nmean: {}\nstd: {}\nmax: {}\n\
                             singletons: {}\nclusters: {}\nseqs: {}\ncompletes: {}')\
                             .format(min_, med, mean, std_dev, max_,
                                     singleton_prop,
                                     features_dict['num_clusters'],
                                     features_dict['num_seqs'],
                                     len(features_dict['completes']))
        plt.ylabel(ylabel_string,
                   rotation='horizontal',
                   horizontalalignment='right')
        
        plt.xlabel('CLUSTER SIZE\n|S| = num species, |C| = num completes')
        
        plt.annotate('*', (len(species), len(features_dict['completes'])),
                     color='r', weight='bold')
        # coordinate annotation offset set to x_range/100 by default
        plt.annotate('(|S|, |C|)',
                     xy=(len(species), len(features_dict['completes'])),
                     xytext=(len(species)+(x_range/100), len(features_dict['completes'])),
                     fontsize=8)
        
        plt.title('MCL Clusters Features: {}'.format(args.mcl_file))

        # SET NUMBER OF HISTOGRAM BINS HERE
        num_bins = 3*len(species) + 1
        
        plt.hist(clstr_size_arr,bins=range(0,num_bins))
        plt.show()

    if args.text:
        print('\ncluster size features\nmin: {}\nmedian: {}\nmean: {}\nstd: {}\nmax: {}\nsingletons: {}\nclusters: {}'
                   .format(min_,med,mean,std_dev,max_,singleton_prop, features_dict['num_clusters']))
        print('\nspecies frequencies in clusters: {}'.format(features_dict['species_frequency']))
        print('\nexcluded species frequencies: {}'.format(features_dict['excluded_species_cts']))
        
    if args.print_bijections:
        print('\nspecies/seqs bijections\n{}'.format(features_dict['species_bijections']))

    if args.print_completes:
        print('\ncompletes\n{}'.format(features_dict['completes']))
        
    # create an MFASTA file from mcl dump file.
    # this will create an mfasta containing
    # only clusters which have all species
    # represented, i.e., "completes"
    if args.fasta_file:
        mfasta = open(args.mcl_file + '.mfasta', 'w')
        fa_dict = fasta.fasta_dict(args.fasta_file)
        for i, complete in enumerate(features_dict['completes']):
            mfasta.write('@cluster' + str(i) + '\n')
            for seq in complete:
                mfasta.write('>' + seq + '\n')
                mfasta.write(fa_dict[seq] + '\n')
        mfasta.write('@END\n')
        mfasta.close()

        
if __name__ == "__main__":
    main()
