from Photon_Sim import Simulation
import matplotlib.pyplot as plt
import numpy as np
import math
import random

sim = Simulation(2.0, 30.0, 3.0, 2.0, 30.0, 3.0, 1.58, 1.0, 1.55, math.pi/4, math.pi/2, detector=2)

def efficiency_histogram(n=1000):
    results = np.zeros(n)
    for i in range (n):
        results[i] = sim.random_test()
        print(i)
    plt.hist(results, bins=20)
    plt.title('Rate of Detection for Randomly Generated Photons')
    plt.xlabel('Fraction of Photons Detected')
    plt.ylabel('Instances')
    print(f'Median: {np.median(results)}')
    print(f'Mean: {np.mean(results)}')
    plt.show()

efficiency_histogram()
