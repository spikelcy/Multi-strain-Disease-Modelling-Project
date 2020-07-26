import math

import networkx as nx
import matplotlib.pyplot as plt
import random
from random import choice
from utils import *
import copy


"""
Multistrain Disease Modelling Project

Author: Spike Lee
Project Supervisor: Dr Nic Geard & Dr Rebecca Chisholm , University of Melbourne.

Description: this file includes the functions used to run the model. Based on work done by previous student and Model
             found in Buckee et al.(2004)
             
Reference: Buckee, C., Danon, L., & Gupta, S. (2007). Host community structure and the maintenance of
pathogen diversity. Proceedings of the Royal Society B: Biological Sciences, 274(1619), 1715-1721
"""

# Check for loss of immunity
def check_immune(node,n,sigma):
    new = node.copy()
    for strain in new.copy():
        p_immune = random.random()
        # print("Rolls %s for loss of immunity" % (p_immune))
        if p_immune < sigma:
            new.remove(strain)
            # print("Node %s loses its immunity to strain %s" % (n,strain))
    return new

# Attempts to infect neighbour
def attempt_infection(s,node,j,G,N,gamma,beta,tao,R):
    # use a list to store the current infections of j
    infections_j = list(G.node[j]['current_infected'])
    infections_j_strains = []
    for strain in infections_j:
        infections_j_strains.append(strain[1])
    # use a list to store the current immune state of j
    immune_state_j = list(G.node[j]['current_immune'])
    immune_j_strains = []
    for strain in immune_state_j:
        immune_j_strains.append(strain[1])
    # If the host has not been infected with the strain and has no cross-immunity to it, infection can happen
    if s not in infections_j and s not in immune_state_j:
        overlap = loci_overlap(s[1], immune_j_strains, N)

        # check if infection takes place
        p = random.random()
        # calculate new probability based on overlap
        vulnerability = (1 - overlap ** (1 / gamma)) ** gamma

        p_infect = beta * vulnerability

        if p < p_infect:
            IMR(s,j,node,G,tao,R)

    return

# Once infection is confirmed, check mutation or recombination or just infection
def IMR(s,j,node,G,tao,R):

    #attempt mutation
    new_strain = mutate(s[1],tao)
    #attemp recombination

    if G.node[j]['current_infected']:
        infected_strain = G.node[j]['current_infected'].copy().pop()

        new_strain = recombine(new_strain,infected_strain[1],R)

    infecting_strain = (s[0],new_strain)
    G.node[j]['newly_infected'].add(infecting_strain)
    return


def attempt_recovery(node,n,mu):
    recover = set()
    for strain in node:
        p = random.random()

        if (p < mu):
            recover.add(strain)

    return recover

def record_total_infected(network, host_infected,num_people):
    for infected_record in host_infected.values():
        infected_record.append(0)
    for n, v in network.nodes.items():
        for strain in v['current_infected']:
            host_infected[strain][-1] += 1
    for infected_record in host_infected.values():
        infected_record[-1] = infected_record[-1] / num_people

def record_total_immune(network, host_immune,num_people):
    for immune_record in host_immune.values():
        immune_record.append(0)
    for n, v in network.nodes.items():
        for strain in v['current_immune']:
            host_immune[strain][-1] += 1
    # for immune_record in host_immune.values():
    #     immune_record[-1] = immune_record[-1] / num_people

def record_com_immune(network, host_immune,num_people,partition):
    # print(host_immune)
    for immune_record in host_immune.values():
        immune_record.append(0)
    for n, v in network.nodes.items():
        if n in partition:
            for strain in v['current_immune']:
                host_immune[strain][-1] += 1
    for immune_record in host_immune.values():
        immune_record[-1] = immune_record[-1] / num_people
    # print(host_immune)

def record_com_infected(network, host_infected,num_people,partition):
    for infected_record in host_infected.values():
        infected_record.append(0)
    for n, v in network.nodes.items():
        if n in partition:
            for strain in v['current_infected']:
                host_infected[strain][-1] += 1
    for infected_record in host_infected.values():
        infected_record[-1] = infected_record[-1] / num_people

def calc_diversity(host_infected):
    entropy = 0
    for freq in host_infected.values():
         if freq[-1] != 0:
            entropy += freq[-1] * math.log(1 / freq[-1])
    diversity = entropy / math.log(len(host_infected))
    return diversity


def hamming(strain1, strain2, n_loci):
    distance = 0
    for i in range(n_loci):
        if strain1[i] != strain2[i]:
            distance += 1
    return distance


def calc_discordance(host_infected, n_loci):
    n = 0
    d = 0
    for strain1, freq1 in host_infected.items():
        # print(strain1)
        for strain2, freq2 in host_infected.items():
            if strain1 != strain2 and freq1[-1] != 0 and freq2[-1] != 0:
                n += hamming(strain1, strain2, n_loci) * freq1[-1] * freq2[-1]
                d += freq1[-1] * freq2[-1]
    try:
        discordance = n / d / n_loci
    except ZeroDivisionError:
        discordance = 0
    return discordance

def calc_div_mean(sum,timesteps):
    return sum/timesteps
def calc_disc_mean(sum,timesteps):
    return sum/timesteps



def model(contacts_per_host, mu, sigma, beta, r, tao, gamma, n_loci, n_nodes,z_out ,n_com,randomness, n_seeds, n_steps,community_type):
    # Attributes used for network
    P = randomness
    num_people = n_nodes
    MAX_TIME_STEPS = n_steps # max time steps
    gamma = gamma # degree of cross-immunity
    beta = beta # probability of infection
    mu = mu  # probability of the host losing infection at a time step( used for testing at the moment)
    sigma = sigma  # probability of host losing immunity at a time step
    tao = tao # probability of  mutation
    R = r # probability of recombination
    div_sum = 0
    N = n_loci # Number of loci
    # testing network functions
    n_communities = n_com

    n_out = z_out
    n_in = 10-n_out
    p_out = n_out/(num_people/n_communities)
    p_in =  n_in/(num_people/n_communities)

    partitions = []

    comm_size = int(num_people / n_communities)
    sizes = []
    for i in range(n_communities):
        sizes.append(comm_size)
    ## change community type
    ## 1 == random
    ## 0 == linear
    if community_type == 1:
        G = nx.random_partition_graph(sizes, p_in, p_out)
    else:
        G = linear_community_network(n_communities,num_people,n_in,n_out)
    for i in range(n_communities):
        G1 = G.graph['partition'][i]
        partitions.append(G1)


    # adding attributes to nodes
    for i in range(n_communities):
        G1 = partitions[i]
        for n, v in G.nodes.items():
            v['current_infected'] = set()
            v['current_immune'] = set()
            v['newly_infected'] = set()
            v['newly_recovered'] = set()
            v['next_immune'] = set()
            if n in G1:
                v['index'] = i


    strain_space = generate_ss(N);
    communities_strain_space = []
    for i in range(n_communities):
        temp =[]
        for strain in strain_space:
            temp2 = []
            temp2.append(i)
            temp2.append(strain)
            temp.append(tuple(temp2))
        communities_strain_space.append(temp)

    # Seed infection
    nodes_per_community = int(n_seeds/n_communities)


    # RANDOM SEEDING
    for i in range(n_communities):

            G1 = partitions[i]
            seed_nodes = random.sample(G1, nodes_per_community)
            for node in seed_nodes:
                strain = choice(communities_strain_space[i])
                # print("Community %s has  strain %s"%(i,strain))
                # print(strain)
                G.node[node]['current_infected'].add(strain)

    # Record per community discordance and diversity if needed
    communities_diversity = []
    # communities_discordance = []
    for i in range(n_communities):
        communities_diversity.append(0)
        # communities_discordance.append(0)

    #Start blank record of host immunity per community
    host_immune_com = []
    for i in range(n_communities):
        host_immune = {}
        for n in range(len(communities_strain_space)):
            for strain in communities_strain_space[n]:
                host_immune[strain] = []
        host_immune_com.append(host_immune)

    # Start blank record of host infected per community
    host_infected_com = []
    for i in range(n_communities):
        host_infected = {}
        for n in range(len(communities_strain_space)):
            for strain in communities_strain_space[n]:
                host_infected[strain] = []
        host_infected_com.append(host_infected)

    # Start blank record of host infected total in network
    total_infected = {}
    for n in range(len(communities_strain_space)):
        for strain in communities_strain_space[n]:
            total_infected[strain] = []

    # Start blank record of host infected total in network
    total_immune = {}
    for n in range(len(communities_strain_space)):
        for strain in communities_strain_space[n]:
            total_immune[strain] = []

    # print(total_immune)
    record_total_infected(G,total_infected,num_people)
    record_total_immune(G,total_immune,num_people)


    for i in range(n_communities):
        # print(i)
        G1 = partitions[i]
        host_immune_temp = copy.deepcopy(host_immune_com[i])
        host_infected_temp = copy.deepcopy(host_infected_com[i])
        record_com_immune(G,host_immune_temp,num_people,G1)
        record_com_infected(G,host_infected_temp,num_people,G1)
        host_immune_com[i] = host_immune_temp
        host_infected_com[i] = host_infected_temp


    for t in range(1,MAX_TIME_STEPS):

        for n, v in G.nodes.items():


        # #check immune status of nodes
        # nodes_immune = [n for n, v in G.nodes(data=True) if v['current_immune'] != []]
        # print("nodes_immune: %s" % nodes_immune)
            if v['current_immune']:
                # print("Node %s immune list is %s"%(n,v['current_immune']))
                # check immunity lost
                v['next_immune'] = check_immune(v['current_immune'],n,sigma)


        # #get a list of infected nodes
            if v['current_infected']:
                # print(n)
                # print(v['current_infected'])
                # if G.node[node]['current_infected'] == []:
                #     print("node %s is empty" % node)
                #     break
                #Attempt to infect neighbours
                for j in G[n]:
                    # print("Node %s is connected to Node %s " % (n, j))
                    infections = v['current_infected'].copy()

                    s = infections.pop()


                    attempt_infection(s,n,j,G,N,gamma,beta,tao,R)
                    # print(j)
                    # print(G.node[j]['next_infected'])
                # Attempts to recover from infection
                # check recovery
                v['newly_recovered'] = attempt_recovery(v['current_infected'],n,mu)



        # update the network attributes
        # print("update all values to next values")
        for n, v in G.nodes.items():
            # Print status of each node attribute before the update
            # print("PRE-UPDATE--")
            # print(n)
            # print('newly_infected: %s' % list(v['newly_infected']))
            # print('newly_recovered:%s' % list(v['newly_recovered']))
            # print('current_infected:%s' % list(v['current_infected']))
            # print('next_immune:%s' % list(v['next_immune']))
            # print('current_immune:%s' % list(v['current_immune']))

            # apply changes
            v['current_immune'] = v['next_immune'] | v['newly_recovered']
            v['current_infected'] = v['current_infected'] - v['newly_recovered']
            v['current_infected'] = v['current_infected'] | v['newly_infected']

            # reset temporary hold fields
            v['next_immune'] = set()
            v['newly_recovered'] = set()
            v['newly_infected'] = set()

            # Print status of each node attribute after the update
            # print("POST-UPDATE--")
            # print(n)
            # print('newly_infected: %s' % list(v['newly_infected']))
            # print('newly_recovered:%s' % list(v['newly_recovered']))
            # print('current_infected:%s' % list(v['current_infected']))
            # print('next_immune:%s' % list(v['next_immune']))
            # print('current_immune:%s' % list(v['current_immune']))

        record_total_infected(G,total_infected,num_people)
        record_total_immune(G,total_immune,num_people)
        for i in range(n_communities):
            G1 = partitions[i]
            host_immune_temp = copy.deepcopy(host_immune_com[i])
            host_infected_temp = copy.deepcopy(host_infected_com[i])
            record_com_immune(G, host_immune_temp, num_people, G1)
            record_com_infected(G, host_infected_temp, num_people, G1)
            host_immune_com[i] = host_immune_temp
            host_infected_com[i] = host_infected_temp


        # Record diversity within communities
        for i in range(n_communities):
            # host_infected = host_infected_com[i]
            diversity = calc_diversity(host_infected_com[i])
            # discordance = calc_discordance(host_infected, n_loci)
            # communities_discordance[i].append(discordance)
            # communities_diversity[i].append(diversity)
            # communities_discordance[i] += discordance
            communities_diversity[i] += diversity


        diversity = calc_diversity(total_infected)
        div_sum += diversity

    mean_div = calc_div_mean(div_sum,MAX_TIME_STEPS)
    for i in range(len(communities_diversity)):
        div = communities_diversity[i]
        mean_div_c = calc_div_mean(div,MAX_TIME_STEPS)
        communities_diversity[i] = mean_div_c

    return host_immune_com, host_infected_com,total_immune, communities_diversity,mean_div



