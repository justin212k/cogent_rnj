#!/usr/bin/env python
"""Relaxed Neighbour joining phylogenetic tree estimation.

Note that by default,
negative branch lengths are reset to 0.0 during the calculations.

This code is primarily based off of Gavin Huttley's nj.py code version 1.1,
however, the algorithm does only a local search for neighbors to join, thus
reducing (theoretically at least) computational time to typically O(N^2 log N)

See, for example:
Relaxed Neighbor Joining: A Fast Distance-Based Phylogenetic Tree
Construction Method, by Jason Evans, Luke Sheneman, James Foster

If this algorithm is the bottleneck of an expensive (in computation time) task,
it may be worthwhile to use clearcut, which is written in c.

"""

import numpy
from random import shuffle, seed
from cogent.core.tree import TreeBuilder
from cogent.phylo.util import distanceDictTo2D


__author__ = "Justin Kuczynski"
__copyright__ = ""
__credits__ = ["Justin Kuczynski", "Gavin Huttley", "Peter Maxwell"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

def rnj(dists, no_negatives=True, randomize=True):
    """Computes a tree using the relaxed neighbor joining method
    
    Arguments:
        - dists: dict of (name1, name2): distance
        - no_negatives: negative branch lengths will be set to 0
        - randomize: the algorithm will search nodes randomly until two
        neighbors are found.
    """
        
    constructor = TreeBuilder(mutable=True).createEdge
    (names, d) = distanceDictTo2D(dists)
    
    nodes = [constructor([], name, {}) for name in names]
    
    while len(nodes) > 2:
        # Eliminate one node per iteration until 2 left
        num_nodes = len(nodes)
        
        # compute r (normalized), the sum of all pairwise distances
        # the normalization is over (num - 2), since later for a given i, j
        # distance(i, j) will be removed, and distance(i, i) = 0 always
        r = numpy.sum(d, 0) * 1./(num_nodes-2.)
        
        # find two nodes i, j that are minimize each other's 
        # transformed distance
        node_indices = range(num_nodes)
        if randomize == True:
            shuffle(node_indices)
        chose_pair = False
        
        # coefficient used calculating transformed distances
        coef = num_nodes * 1./(num_nodes - 2.)
        for i in node_indices:
        # find i's closest, call it j
        
            # xformed_dists is a list of T_i,j for all j
            xformed_dists = coef*d[i] - r - r[i]
        
            # give distance to self a bogus but nonminimum value
            xformed_dists[i] = numpy.abs(xformed_dists[0])*2. +\
                numpy.abs(xformed_dists[num_nodes - 1])*2.
        
            j = numpy.argmin(xformed_dists)
        
        
        # now find j's closest
            xformed_dists = coef*d[j] - r - r[j]
            xformed_dists[j] = numpy.abs(xformed_dists[0])*2. +\
                numpy.abs(xformed_dists[num_nodes - 1])*2.
            
            # if i and j are each other's minimum, choose this (i, j) pair
            if i == numpy.argmin(xformed_dists):
                # choose these i, j
                chose_pair = True
                break
        
        if not chose_pair:
            raise Exception("didn't choose a pair of nodes correctly")
        assert i != j, (i, j)
        
        # Branch lengths from i and j to new node
        nodes[i].Length = 0.5 * (d[i,j] + r[i] - r[j])
        nodes[j].Length = 0.5 * (d[i,j] + r[j] - r[i])
            
        # no negative branch lengths
        if no_negatives:
            nodes[i].Length = max(0.0, nodes[i].Length)
            nodes[j].Length = max(0.0, nodes[j].Length)
        
        # Join i and k to make new node
        new_node = constructor([nodes[i], nodes[j]], None, {})
        
        # Store new node at i
        new_dists = 0.5 * (d[i] + d[j] - d[i,j])
        d[:, i] = new_dists
        d[i, :] = new_dists
        d[i, i] = 0.0
        nodes[i] = new_node
        
        # Eliminate j
        d[j, :] = d[num_nodes-1, :]
        d[:, j] = d[:, num_nodes-1]
        assert d[j, j] == 0.0, d
        d = d[0:num_nodes-1, 0:num_nodes-1]
        nodes[j] = nodes[num_nodes-1]
        nodes.pop()
    
    # no negative branch lengths
    if len(nodes[0].Children) < len(nodes[1].Children):
        nodes.reverse()
    
    # 2 remaining nodes will be [root, extra_child]
    nodes[1].Length = d[0,1]
    if no_negatives:
        nodes[1].Length = max(0.0, nodes[1].Length)
    
    #Need to replace nodes[0] with new root
    nodes[1].Parent = nodes[0]
    return constructor(nodes[0].Children, 'root', {}).deepcopy()
