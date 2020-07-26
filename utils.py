import random
import itertools
from random import sample, choice
import networkx as nx

"""
Utilities

Author: Spike Lee
Project Supervisor: Dr Nic Geard & Dr Rebecca Chisholm , University of Melbourne.

Description: this file contains general utilises used by the other files in this project.    
"""

#Generate all possible strains based on number of loci(N).
def generate_ss(N):
    return ["".join(seq) for seq in itertools.product("01", repeat=N)]

#Check if current strain has overlap in immune memory.
def loci_overlap(strain, immune_memory, N):
    i = 0  # keep track of index for loci
    num_of_same_loci = 0
    while i < len(strain):
        # traverse the immune list to find same loci at index i
        for s in immune_memory:
            if strain[i] == s[i]:
                num_of_same_loci += 1
                break
        i += 1
    # Might not use fraction
    overlap = num_of_same_loci / N

    return overlap

# Mutates input strain based on probability(p)
def mutate(strain,p):
    new_strain = ''
    for bit in strain:
        if random.random() < p:
            #use mod to represent bits
            change_bit = (int(bit) + 1) % 2
            new_strain += str(change_bit)
        else:
            new_strain += bit
    # print("%s mutates to %s" % (strain, new_strain))
    return new_strain

# Take two strains to produce a child strain
def recombine(strain_1,strain_2,p):
    new_strain = ''
    #iterate through length of strains
    for i in range(len(strain_1)):

        # check if it will combine
        if random.random() < p:
            # Recombine by choosing bit
            change_bit = choice([strain_1[i],strain_2[i]])
            # print("%s and %s combines to form %s" %(strain_1[i], strain_2[i], str(change_bit)))
            new_strain += str(change_bit)
        # if no recombination at allele, keep bit from infecting strain ( follows peilin)
        else:
            new_strain += strain_1[i]
    return new_strain

def find_discordant_strain(strain):
    discordant_strain = ''
    n = 0
    while n <len(strain):
	    discordant_strain += str(1 - int(strain[n]))
	    n += 1
    return discordant_strain


# Create Linear community network
def linear_community_network(n_communities,num_people,z_in,z_out):
    comm_size = int(num_people / n_communities)
    sizes = []
    n_out = z_out
    n_in = z_in
    p_out = n_out / (num_people / n_communities)
    p_in = n_in / (num_people / n_communities)
    for i in range(n_communities):
        sizes.append(comm_size)
    nodelist = range(0, sum(sizes))
    partitions = []
    print(n_in)
    print(n_out)
    g = nx.Graph()
    size_cumsum = [sum(sizes[0:x]) for x in range(0, len(sizes) + 1)]
    g.graph['partition'] = [set(nodelist[size_cumsum[x]:size_cumsum[x + 1]])
                            for x in range(0, len(size_cumsum) - 1)]

    # Setup nodes and graph name
    for block_id, nodes in enumerate(g.graph['partition']):
        for node in nodes:
            g.add_node(node, block=block_id)

    for i in range(n_communities):
        G1 = g.graph['partition'][i]
        partitions.append(G1)

    for i in range(n_communities):
        G1 = partitions[i]
        if i == 0:
            for node in G1:
                for node2 in G1:
                    if node != node2:
                        p = random.random()
                        if p < p_in:
                            g.add_edge(node, node2)
                for noderight in partitions[i + 1]:
                    p = random.random()
                    if p < p_out:
                        g.add_edge(node, noderight)
        elif (i == n_communities - 1):
            for node in G1:
                for node2 in G1:
                    if node != node2:
                        p = random.random()
                        if p < p_in:
                            g.add_edge(node, node2)
            for node in G1:
                for nodeleft in partitions[i - 1]:
                    p = random.random()
                    if p < p_out:
                        g.add_edge(node, nodeleft)
        else:
            for node in G1:
                for node2 in G1:
                    if node != node2:
                        p = random.random()
                        if p < p_in:
                            g.add_edge(node, node2)
            for node in G1:
                for nodeleft in partitions[i - 1]:
                    p = random.random()
                    if p < p_out:
                        g.add_edge(node, nodeleft)
                for noderight in partitions[i + 1]:
                    p = random.random()
                    if p < p_out:
                        g.add_edge(node, noderight)
    return g

