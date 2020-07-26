import json
from multiprocessing.pool import Pool


import numpy as np
from smt.sampling_methods import LHS
import time

from network import model

"""
Latin HyperCube Sampling(LHS) of Models

Author: Spike Lee
Project Supervisor: Dr Nic Geard & Dr Rebecca Chisholm , University of Melbourne.

Description: this file uses LHS to create model parameters to run with the created multi-strain disease models. Parameter
             range can be changed and results are saved in text files for plotting.     
"""

if __name__ == '__main__':
    start_time = time.time()
    # parameter range (changed to same format)
    CONTACTS_PER_HOST = [10, 10]
    MU = [1/7, 1/7]  # recovery probability
    SIGMA = [1/23, 1/23]  # immunity lost probability
    BETA = [0.2472, 0.2472]  # infection probability
    R = [0.0953, 0.0953]  # recombination probability per allele
    TAO = [0.0042, 0.0042]  # mutation probability per allele
    GAMMA = [4, 4]  # cross-immunity
    N_LOCI = [2,2]
    N_NODES = [256, 256]
    Z_OUT = [0.5,0.5]
    N_COM = [3,10]

    # fixed parameters
    LINEAR = 0  # host contact network structure
    RANDOM = 1
    N_SEEDS = 12
    N_STEPS = 1000

    # create multi-dimension parameter space
    parameter_space = np.array([
        CONTACTS_PER_HOST,
        MU,
        SIGMA,
        BETA,
        R,
        TAO,
        GAMMA,
        N_LOCI,
        N_NODES,
        Z_OUT,
        N_COM
    ])

    # LHS
    sampling = LHS(xlimits=parameter_space)
    n_sample_points = 160
    parameter_samples = sampling(n_sample_points)

    parameter_sample_list = parameter_samples.tolist()
    # used to save div and disc for scatter plot
    ran_div = []
    lin_div = []
    data = {'ran_data': []}
    data_linear = {'linear_data': []}
    pool = Pool(processes=2)

    # round int parameters
    print("Running...")
    for sample in parameter_sample_list:
        sample[0] = int(round(sample[0]))
        sample[7] = int(round(sample[7]))
        sample[10] = int(round(sample[10]))
        sample[8] = sample[10]*100
        print(sample)
        print(parameter_sample_list.index(sample))
        host_immune_fixed = []
        host_immune2_fixed = []
        results_linear = pool.apply_async(model, (*sample, 1, N_SEEDS, N_STEPS, LINEAR))
        results_random = pool.apply_async(model, (*sample, 1, N_SEEDS, N_STEPS, RANDOM))

        host_immune, host_infected, total_immune, com_diversity, total_diversity = results_linear.get()
        host_immune2, host_infected2, total_immune2, com_diversity2, total_diversity2 = results_random.get()


        data_linear['linear_data'].append({
            'CONTACTS_PER_HOST': sample[0],
            'MU': sample[1],
            'SIGMA': sample[2],
            'BETA': sample[3],
            'R': sample[4],
            'TAO': sample[5],
            'GAMMA': sample[6],
            'N_LOCI': sample[7],
            'N_NODES': sample[8],
            'Z_OUT':sample[9],
            'N_COM':sample[10],
            'com_diversity': com_diversity,
            'total_diversity': total_diversity
        })
        data['ran_data'].append({
            'CONTACTS_PER_HOST': sample[0],
            'MU': sample[1],
            'SIGMA': sample[2],
            'BETA': sample[3],
            'R': sample[4],
            'TAO': sample[5],
            'GAMMA': sample[6],
            'N_LOCI': sample[7],
            'N_NODES': sample[8],
            'Z_OUT': sample[9],
            'N_COM': sample[10],
            'com_diversity': com_diversity2,
            'total_diversity': total_diversity2
        })

    with open('results_linear_LHS.txt', 'w') as outfile:
        json.dump(data_linear, outfile)
    with open('results_random_LHS.txt', 'w') as outfile:
        json.dump(data, outfile)

    print("--- %s seconds ---" % (time.time() - start_time))
