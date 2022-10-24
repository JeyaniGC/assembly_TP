#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import textwrap

__author__ = "Jeyani George Clement"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Jeyani George Clement"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Jeyani George Clement"
__email__ = "jeyanigeorgeclement@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
       
       Parameters:
        path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Read dasta file to create a sequence geeerator.
    
    Parameters
    ----------
    fastq_file : Fastq file containing the sequence it's id and it's quality.
    
    Returns
    -------
    A sequences generator.
    """
    with open(fastq_file, "r") as filin:
        for i in filin:
            yield next(filin).strip()
            next(filin)
            next(filin)


def cut_kmer(read, kmer_size):
    """Read a sequence to create a k-mer generator.
    
    Parameters
    ----------
    read : str()
        Represents the sequence to be cut.
    kmer_size : int()
        Size of the cut in the sequence.

    Returns
    -------
    A k-mer generator.
    """
    
    for i in range(len(read)-(kmer_size-1)):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Creation of a dictionary whose keys are the different k-mers and whose 
    value is the occurrence of these k-mers.
    
    Parameters
    ----------
    fastq_file : A fastq file
        Contains the id of a sequence, the sequence and its quality.
    kmer_size : int()
        Size of the cut in the sequence.

    Returns
    -------
    A dictionary listing all k-mers and their occurrences within the fastq file.
    """
    seq = "".join(list(read_fastq(fastq_file)))
    kmer = cut_kmer(seq, kmer_size)
    dict_kmer = {}
    
    for node in kmer:
        if node in dict_kmer:
            dict_kmer[node] += 1
        else:
            dict_kmer[node] = 1
            
    return dict_kmer


def build_graph(kmer_dict):
    
    
    g = nx.DiGraph()
    
    for node, weight in kmer_dict.items():
        g.add_edge(node[:-1], node[1:], weight =  weight)
        
    return g


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node and delete_sink_node == False:
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node and delete_entry_node == False:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    
    return graph

def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    
    weight_list = std(weight_avg_list)

    if weight_list > 0:
        best_path_index = weight_avg_list.index(max(weight_avg_list))
    elif weight_list == 0:
        length_std = std(path_length)
        if length_std > 0:
            best_path_index = path_length.index(max(path_length))
        elif length_std == 0:
            best_path_index = randint(0,len(path_list)-1)

    for path in path_list:
        if path != path_list[best_path_index]:
            graph = remove_paths(graph, [path], delete_entry_node, delete_sink_node)
    
    return graph

def path_average_weight(graph, path):
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    nodes = list(graph.nodes())
    starting_nodes = []
    
    for node in nodes: 
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
            
    return starting_nodes

def get_sink_nodes(graph):
    nodes = list(graph.nodes())
    output_nodes = []
    
    for node in nodes: 
        if len(list(graph.successors(node))) == 0:
            output_nodes.append(node)
            
    return output_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    
    for start in starting_nodes:
        for end in ending_nodes: 
            if nx.has_path(graph, start, end):
                for simple_path in nx.all_simple_paths(graph, start, end):
                    seq_contig = simple_path[0]
                    for node in simple_path[1:]:
                        seq_contig += node[-1]
                    contigs.append(tuple((seq_contig, len(seq_contig))))
    return contigs

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as filout:
        for i in range(len(contigs_list)):
            filout.write(f">contig_{i} len={contigs_list[i][1]}\n")
            filout.write(f"{textwrap.fill(contigs_list[i][0], width = 80)}\n")   
            
def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
