from Photon_Sim import Simulation
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math
import random

sim = Simulation(2.0, 30.0, 3.0, 2.0, 30.0, 2.0, 1.58, 1.0, 1.55, detector=2)

def random_efficiency_histogram(n=1000):
    # Creates detection rate histogram for photons randomly generated inside the scintillator
    # Currently no support for light pipe or airgap case

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

def efficiency_histogram(*args, n=1000, plot=True):
    # Creates detection rate histogram for n=1000 photons generated along a random electron intersection path
    # inside the scintillator, allowing for modular addition of light pipes and other mediums
    # Args = (2d array: dimensions of every intermediate medium between the scintillator and detector
    #                   shape = (number of mediums, 3 (dimensions l, w, h))
    #         n floats: index of refraction for each intermediate medium)
    # Returns = (float: median of results
    #            array: detection rates)

    results = np.zeros(n)
    y = np.random.uniform(-1.0, 1.0, n) * sim.w / 2
    z = np.random.uniform(-1.0, 1.0, n) * sim.h / 2
    for i in range (n):
        results[i] = sim.run(y[i], z[i], *args)
        print(i)
    median = np.median(results)
    #print(f'Median: {median}')
    #print(f'Mean: {np.mean(results)}')
    if plot:
        fig = plt.figure()
        plt.hist(results, bins=20)
        plt.title('Rate of Detection for Photons from\nRandomly Generated Electron Paths')
        plt.xlabel('Fraction of Photons Detected'   )
        plt.ylabel('Instances')
        plt.show(block=False)
        plt.pause(.1)
    return median, results

def gap_efficiency_scatter(point=101):
    # Scatterplot of detection rate by length of the airgap between the light and the SiPM window
    # Point number of points are from 0 to 1 mm and one tenth of that are from 1 to 2 mm
    # Best fit is currently linear, though that is visibly not the best model

    ds = np.concatenate((np.linspace(0, 1, point), np.linspace(1.1, 2, point // 10)))
    eff = np.zeros_like(ds)
    err = np.zeros_like(ds)
    for i in range(len(ds)):
        dimensions = np.array([[2.0, 0.125, 3.0], [2.0, 54.86, 3.0], [100.0, ds[i], 100.0]])
        eff[i], results = efficiency_histogram(dimensions, 1.57, 1.502, 1.0, plot=False)
        err[i] = np.std(results)
        print(f'{round(100 * i / (point + point // 10), 1)}%')
    b, m = np.polynomial.polynomial.Polynomial.fit(ds, eff, 1).convert().coef
    fig = plt.figure()
    plt.errorbar(ds, eff, yerr=err, fmt='o', ls='', elinewidth=.5)
    plt.plot(ds, (m * ds) + b, label=f'Best Fit: \n Slope: {m:.2f} \n Intercept: {b:.2f}')
    plt.title('Rate of Detection by Gap Distance between Light Pipe and SiPM Window')
    plt.xlabel('Airgap Distance (mm)')
    plt.ylabel('Rate of Photon Detection')
    plt.legend()
    plt.show(block=False)
    plt.pause(1)

def pipe_length_efficiency_scatter(point=101):
    # Scatterplot of detection rate by length of the light pipe
    # Point number of points are from 50 to 60 mm
    # Best fit is currently linear, though that is likely not the best model

    ds = np.linspace(50, 60, point)
    eff = np.zeros_like(ds)
    err = np.zeros_like(ds)
    for i in range(len(ds)):
        dimensions = np.array([[2.0, 0.125, 3.0], [2.0, ds[i], 3.0]])
        eff[i], results = efficiency_histogram(dimensions, 1.57, 1.502, plot=False)
        err[i] = np.std(results)
        print(f'{round(100 * i / point, 1)}%')
    b, m = np.polynomial.polynomial.Polynomial.fit(ds, eff, 1).convert().coef
    fig = plt.figure()
    plt.errorbar(ds, eff, yerr=err, fmt='o', ls='', elinewidth=.5)
    plt.plot(ds, (m * ds) + b, label=f'Best Fit: \n Slope: {m:.2f} \n Intercept: {b:.2f}')
    plt.title('Rate of Detection by Light Pipe Length')
    plt.xlabel('Light Pipe Length (mm)')
    plt.ylabel('Rate of Photon Detection')
    plt.legend()
    plt.show(block=False)
    plt.pause(1)

def absorption_histogram(n=1000):
    # Creates histograms of the length travelled by a photon and recursive depth before it escapes,
    # attenuates, or is detected.

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

def heat_map(*args, run=sim.run):
    # Heatmap of detection from photons generated along an electron beams intersection line with
    # the scintillator at different points
    # Args = (2d array: dimensions of every intermediate medium between the scintillator and detector
    #                   shape = (number of mediums, 3 (dimensions l, w, h))
    #         n floats: index of refraction for each intermediate medium)

    y = np.linspace(-sim.w / 2, sim.w / 2, 100)
    z = np.linspace(-sim.h / 2, sim.h / 2, 10)
    heat = np.zeros((z.size, y.size))
    for i in range(len(y)):
        for j in range(len(z)):
            #print('{i} {j}')
            heat[j, i] = run(y[i], z[j], *args)
        print(i)

    fig, ax = plt.subplots()

    im = ax.imshow(heat, cmap='magma', extent=(-sim.w / 2, sim.w / 2, -sim.h / 2, sim.h / 2))
    plt.colorbar(im, ax=ax, label="Rate of Detection")
    ax.set_title('Detection Rate Heatmap')
    ax.invert_xaxis()
    plt.show(block=False)
    plt.pause(.1)

def paths_display(*args, sample=100,
                  dimensions=np.array([[2.0, 0.125, 3.0], [2.0, 54.86, 3.0], [100.0, 0.1, 100.0]])):
    # Displays all paths and end points of photons generated from an intersection with the center of the scintillator
    # Currently, changing which cases are displayed requires manual tweaking in the run method
    # Args = (n floats: index of refraction for each intermediate medium)

    sim.history = True
    print(sim.run(0, 0, dimensions, *args, n=sample))
    fig, ax = plt.subplots(figsize=(30, 5))
    # print(f'First TIR: {math.asin(sim.n2 / sim.n1)}\nSecond TIR: {math.asin(sim.n2 / args[1])}')

    for ls in sim.paths:
        path = np.array(ls, dtype=float)
        y_offset = 0
        y_contin = []
        for i in range(len(path)):
            if i > 0 and np.isclose(path[i, 0], path[i - 1, 0]) and \
                    np.isclose(path[i, 2], path[i - 1, 2]) and not np.isclose(path[i, 1], path[i - 1, 1]):
                # segment boundary!
                y_offset += (path[i - 1, 1] - path[i, 1])
            y_contin.append([path[i, 0], path[i, 1] + y_offset, path[i, 2]])
        y_contin = np.array(y_contin)
        plt.plot(y_contin[:, 1], y_contin[:, 2], alpha=0.4, lw=0.5)
        plt.scatter(y_contin[-1, 1], y_contin[-1, 2], zorder=10000,
                    color=('black' if abs(y_contin[-1,1] - (sim.w/2 + dimensions[0, 1] + dimensions[1, 1])) < 1e-2 else
                           'green' if abs(y_contin[-1,2]) == float(dimensions[1, 2])/2 else
                           'blue' if abs(y_contin[-1, 0]) == float(dimensions[1, 0]) / 2
                           else 'red')
                    )
        #print(y_contin)

    # path = np.array(sim.paths[2])
    # for i in range(path.shape[0] - 1):
    #     if path[i, 0] == path[i + 1, 0]:
    #         path[i + 1:, 1] += (path[i, 1] - path[i + 1, 1])
    # print(path)
    # plt.plot(path[:, 1], path[:, 2], label=f'Successful Paths')
    plt.title('Path of a Detected Photon in Axes Normal to Electron Axis (x)')
    plt.xlabel('y')
    plt.ylabel('z')
    plt.ylim(-2, 2)
    plt.xlim(-16.5, 72)
    scin = patches.Rectangle((-sim.w / 2, -sim.h / 2), sim.w, sim.h,
                             linewidth=1, edgecolor='cyan', facecolor='cyan', label='Scintillator', zorder=0)
    pipe = patches.Rectangle(((sim.w / 2) + float(dimensions[0, 1]), -float(dimensions[1, 2]) / 2),
                             float(dimensions[1, 1]), float(dimensions[1, 2]),
                             linewidth=1, edgecolor='xkcd:lemon', facecolor='xkcd:lemon', label='Light Pipe', zorder=0)
    ax.add_patch(scin)
    ax.add_patch(pipe)
    plt.legend()
    plt.show(block=False)
    plt.pause(1)


# Efficiency histogram and statistics for no pipe, 2x2 pipe, and 3x2 pipe scenarios:
# dimensions = np.array([[100.0, 0.1, 100.0]])
# median, results = efficiency_histogram(dimensions, 1.0)
# print(f'Median: {median}')
# print(f'Mean: {np.mean(results)}')
# print(f'Standard Deviation: {np.std(results)}')
#
# dimensions = np.array([[2.0, 0.125, 2.0], [2.0, 54.86, 2.0], [100.0, 0.1, 100.0]])
# median, results = efficiency_histogram(dimensions, 1.57, 1.502, 1.0)
# print(f'Median: {median}')
# print(f'Mean: {np.mean(results)}')
# print(f'Standard Deviation: {np.std(results)}')
#
# dimensions = np.array([[2.0, 0.125, 3.0], [2.0, 54.86, 3.0], [100.0, 0.1, 100.0]])
# median, results = efficiency_histogram(dimensions, 1.57, 1.502, 1.0)
# print(f'Median: {median}')
# print(f'Mean: {np.mean(results)}')
# print(f'Standard Deviation: {np.std(results)}')

#dimensions = np.array([[2.0, 0.125, 3.0], [2.0, 54.86, 3.0]])
#efficiency_histogram(dimensions, 1.57, 1.502)
# absorption_histogram(10000)



# Heat map for no pipe, 3x2 pipe, and 2x2 pipe scenarios:
# dimensions = np.array([[100.0, 0.1, 100.0]])
# heat_map(dimensions, 1.0, run=sim.run)
# dimensions = np.array([[2.0, 0.125, 3.0], [2.0, 54.86, 3.0], [100.0, 0.1, 100.0]])
# heat_map(dimensions, 1.57, 1.502, 1.0, run=sim.run)
# dimensions = np.array([[2.0, 0.125, 2.0], [2.0, 54.86, 2.0], [100.0, 0.1, 100.0]])
# heat_map(dimensions, 1.57, 1.502, 1.0, run=sim.run)

#gap_efficiency_scatter()
#pipe_length_efficiency_scatter()

#Displays all paths for a 2x2 light pipe
# dimensions = np.array([[2.0, 0.125, 2.0], [2.0, 54.86, 2.0], [100.0, 0.1, 100.0]])
# paths_display(1.57, 1.502, 1.0, dimensions=dimensions, sample=10000)

plt.show()