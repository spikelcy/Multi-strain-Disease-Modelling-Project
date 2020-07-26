import networkx as nx
import matplotlib.pyplot as plt
import random

"""
Linear communities Network testing and visualization

Author: Spike Lee
Project Supervisor: Dr Nic Geard & Dr Rebecca Chisholm , University of Melbourne.

Description: this file allows a user to create and plot networks with linearly connected communities based on their 
             predefined input values. Definitions and descriptions of attributes can be seen in reference papers and
             final thesis paper.            
"""
n_communities = 4
num_people = 40
mean_contacts = 10
p_out = 0.05 # the network is random if its 0.05
n_out = int((num_people / n_communities) * p_out) + 1
n_in = mean_contacts - n_out
p_in = n_in / (num_people / n_communities)

comm_size = int(num_people/n_communities)
sizes = []
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
                        g.add_edge(node,node2)
            for noderight in partitions[i+1]:
                p = random.random()
                if p < p_out:
                    g.add_edge(node, noderight)
    elif(i == n_communities-1):
        for node in G1:
            for node2 in G1:
                if node != node2:
                    p = random.random()
                    if p < p_in:
                        g.add_edge(node,node2)
        for node in G1:
            for nodeleft in partitions[i-1]:
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


nx.draw(g)
plt.show()