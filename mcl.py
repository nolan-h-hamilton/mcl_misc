"""Compute some basic features of MCL results

    Basic Usage Example:
        # produce a graphic with a log-scaled histogram of cluster sizes and various other features
        python mcl.py data/mcl_sample.out --species "AAMBTRI, ANUPADV, FEQUDIF, ALIRTUL, GPINTAE" --plot

        # produce a graphic with a log-scaled histogram of cluster sizes and various other features
        # normalize flag causes species frequency dict (see cluster_features()) to be normed by the
        # total number of clusters: "what fraction of clusters contained AAMBTRI?"
        python mcl.py data/mcl_sample.out --species "AAMBTRI, ANUPADV, FEQUDIF, ALIRTUL, GPINTAE" --plot --normalize
       
        # print some basic cluster features to terminal
        python mcl.py data/mcl_sample.out --species "AAMBTRI, ANUPADV, FEQUDIF, ALIRTUL, GPINTAE" --text
        
        # normalize, get bijections, print and plot
        python mcl.py 06042020.mci.I14.out.dump --species "AAMBTRI, ANUPADV, ALIRTUL, FEQUDIF, GPINTAE, GCYCMIC, AKADHET, AMUSACU, AGLYMAX, ASOLPEN, APRUPER, FLYGJAP, ASOLPEN, ABETVUL, ABRADIS" --text --normalize --print_bijections --plot
"""

import sys
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt


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
            'species_bijections': a list containing clusters with len(species) sequences with 
                 every species represented--clusters that exhibit a one-to-one mapping between
                 sequences and species
            'excluded_species_cts': a dictionary containing species keys and how often the species'
                 absence from a cluster was responsible for the cluster's disqualification as a "bijection"
             
    Notes:
        This function should provide all data needed from the mcl file
            in a single traversal. If more features are of interest, this
            function should be modified to compute them without additional
            traversals (i.e., computed in the main loop) if possible.
    """
    
    clusters = []
    species_frequency = dict.fromkeys(species, 0)
    excluded_species_cts = dict.fromkeys(species, 0)
    cluster_sizes = []
    singletons = []
    num_clusters = 0
    size_distr_dict = {'mean_clstr_size': 0.0, 'min_clstr_size': 1000.0,
                       'max_clstr_size': 0.0}
    species_bijections = []
    
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
                # note that the species of a sequence is assumed
                # to be present in the first 7 characters of the
                # sequence name
                present_specs = [x[0:7] for x in cluster]
                set_species = set(species)
                set_pres = set(present_specs)
                spec_diff = set_species - set_pres
                
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
            'num_clusters': num_clusters,
            'species_bijections': species_bijections,
            'excluded_species_cts': excluded_species_cts,
            'clusters': clusters}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mcl_file', help='file containing MCL results')
    parser.add_argument('--species', help='comma-separated species names, ex: "ANUPADV, AAMBTRI"',)
    parser.add_argument('--plot', help='display histogram of clusters features',
                        action='store_true')
    parser.add_argument('--text', help='print cluster features to terminal',
                        action='store_true')
    parser.add_argument('--normalize', help='normalize species frequencies',
                        action='store_true')
    parser.add_argument('--print_bijections', help='print a list of clusters that are one-to-one-mappings',
                        action='store_true')
    
    args = parser.parse_args()
    
    features_dict = cluster_features(args.mcl_file, args.species.split(', '),
                                      normalize=args.normalize)

    clstr_size_arr = np.array(features_dict['cluster_sizes'])
    mean = round(features_dict['size_distr_dict']['mean_clstr_size'],2)
    max_ = round(features_dict['size_distr_dict']['max_clstr_size'],2)
    min_ = round(features_dict['size_distr_dict']['min_clstr_size'],2)
    singleton_prop = round(len(features_dict['singletons'])
                           / float(features_dict['num_clusters']),2)
    # one-pass std devation estimators may not work for larger cluster sizes, use traditional
    # method
    std_dev = round(np.std(clstr_size_arr),2)
    med = round(np.median(clstr_size_arr),2)

    if args.plot:
        plt.yscale('log')
        plt.ylabel('min: {}\nmedian: {}\nmean: {}\nstd: {}\nmax: {}\nsingletons: {}\nclusters: {}'
                   .format(min_,med,mean,std_dev,max_,singleton_prop, features_dict['num_clusters']),
                   rotation='horizontal',
                   horizontalalignment='right')
        plt.xlabel(features_dict['species_frequency'])
        plt.title('MCL Clusters Features')
        plt.hist(clstr_size_arr)
        plt.show()

    if args.text:
        print('\ncluster size features\nmin: {}\nmedian: {}\nmean: {}\nstd: {}\nmax: {}\nsingletons: {}\nclusters: {}'
                   .format(min_,med,mean,std_dev,max_,singleton_prop, features_dict['num_clusters']))
        print('\nspecies frequencies in clusters: {}'.format(features_dict['species_frequency']))
        print('\nexcluded species frequencies: {}'.format(features_dict['excluded_species_cts']))
        
    if args.print_bijections:
        print('\nspecies/seqs bijections\n{}'.format(features_dict['species_bijections']))
    
if __name__ == "__main__":
    main()
