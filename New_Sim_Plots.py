from Photon_Sim import Simulation
import matplotlib.pyplot as plt
import numpy as np
import math
import random

sim = Simulation(2.0, 30.0, 3.0, 2.0, 30.0, 2.0, 1.58, 1.0, 1.55, detector=2)

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
    # plt.hist(hits[:, 1], bins=30)
    # plt.title('Number of Ray Trace Recursions Before Detection')
    # plt.xlabel('Number of Reflections')
    # plt.ylabel('Instances')
    # plt.yscale('log')
    # plt.show()
    plt.hist(misses[:, 0], bins=30)
    plt.title('Lengths Photons Traveled Before \n Escape/Absorption/Geometric Impossibility')
    plt.xlabel('length (mm)')
    plt.ylabel('Instances')
    plt.yscale('log')
    plt.show()
    # plt.hist(misses[:, 1], bins=30)
    # plt.title('Number of Ray Trace Recursions Before \n Escape/Absorption/Geometric Impossibility')
    # plt.xlabel('Number of Reflections')
    # plt.ylabel('Instances')
    # plt.yscale('log')
    # plt.show()

def heat_map(*args, run=sim.input_test):
    y = np.linspace(-sim.w / 2, sim.w / 2, 100)
    z = np.linspace(-sim.h / 2, sim.h / 2, 10)
    heat = np.zeros((z.size, y.size))
    for i in range(len(y)):
        for j in range(len(z)):
            #print('{i} {j}')
            heat[j, i] = run(y[i], z[j], *args, n=1000)
        print(i)

    fig, ax = plt.subplots()

    im = ax.imshow(heat, cmap='magma', extent=(-sim.w / 2, sim.w / 2, -sim.h / 2, sim.h / 2))
    plt.colorbar(im, ax=ax, label="Rate of Detection")
    ax.set_title('Detection Rate Heatmap')
    ax.invert_xaxis()
    plt.show()

# efficiency_histogram()

# absorption_histogram(10000)
dimensions = np.array([[2.0, 0.125, 3.0], [2.0, 54.86, 3.0], [100.0, 0.1, 100.0]])
#print(sim.run(0, 0, dimensions, 1.57, 1.502, 1.0))
heat_map(dimensions, 1.57, 1.502, 1.0, run=sim.run)
#heat_map()