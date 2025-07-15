from Photon_Sim import Simulation
import matplotlib.pyplot as plt
import numpy as np
import math
import random

sim = Simulation(2.0, 30.0, 3.0, 2.0, 30.0, 2.0, 1.58, 1.0, 1.55, math.pi/4, math.pi/2, detector=2)

def efficiency_histogram(n=1000):
    results = np.zeros(n)
    for i in range (n):
        results[i] = sim.random_test()[0]
        print(i)
    plt.hist(results, bins=20)
    plt.title('Rate of Detection for Randomly Generated Photons')
    plt.xlabel('Fraction of Photons Detected')
    plt.ylabel('Instances')
    print(f'Median: {np.median(results)}')
    print(f'Mean: {np.mean(results)}')
    plt.show()

def absorption_histogram(n=1000):
    sim.iterations = n
    _, hits, misses = sim.random_test()
    hits = np.array(hits)
    misses = np.array(misses)
    print(hits.shape)
    plt.hist(hits[:, 0], bins=30)
    plt.title('Lengths Photons Traveled Before Detection')
    plt.xlabel('length (mm)')
    plt.ylabel('Instances')
    plt.yscale('log')
    plt.show()
    plt.hist(hits[:, 1], bins=30)
    plt.title('Number of Ray Trace Recursions Before Detection')
    plt.xlabel('Number of Reflections')
    plt.ylabel('Instances')
    plt.yscale('log')
    plt.show()
    plt.hist(misses[:, 0], bins=30)
    plt.title('Lengths Photons Traveled Before \n Escape/Absorption/Geometric Impossibility')
    plt.xlabel('length (mm)')
    plt.ylabel('Instances')
    plt.yscale('log')
    plt.show()
    plt.hist(misses[:, 1], bins=30)
    plt.title('Number of Ray Trace Recursions Before \n Escape/Absorption/Geometric Impossibility')
    plt.xlabel('Number of Reflections')
    plt.ylabel('Instances')
    plt.yscale('log')
    plt.show()

# efficiency_histogram()

absorption_histogram(10000)