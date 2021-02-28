from Photon_Sim import *
import matplotlib.pyplot as plt
import numpy as np
import csv

def n1_vary(wp, hp, air_gap):
    n1 = []
    eff = []

    for i in range(11):
        sim = Simulation(20.0, 3.0, 2.0, 0.0, wp, hp, 1.0 + (i / 10), 1.0, 1, air_gap)
        n1.append(sim.n1)
        eff.append(sim.run())

    fig, ax = plt.subplots()

    ax.plot(n1, eff)


    ax.set(xlabel='Refractive Index of Bar (n1)', ylabel='Light Efficiency (eff)',
           title='Light Efficiency as a Function of Refractive Index')
    ax.set_ylim([0, 100])
    ax.set_xlim([1.0, 2.0])

    plt.text(1.05, 90, "l = 20.0")
    plt.text(1.05, 85, "w = 3.0")
    plt.text(1.05, 80, "h = 2.0")
    plt.text(1.05, 75, "wp = " + str(wp))
    plt.text(1.05, 70, "hp = " + str(hp))
    plt.text(1.05, 65, "Detector = 1")
    plt.text(1.05, 60, "Airgap = " + str(air_gap))

    plt.show()
    plt.savefig("n1_vary.png")


def tc_vary(wp, hp, air_gap):
    tc = []
    eff = []

    for i in range(25):
        sim = Simulation(50.0, 3.0, 2.0, 0.0, wp, hp, 1.0 + (i / 10), 1.0, 1, air_gap)
        tc.append(math.asin(1/sim.n1)*(180/math.pi))
        eff.append(sim.run())

    fig, ax = plt.subplots()

    if air_gap == False:
        tc_1 = np.arange(0, 46)
        tc_2 = np.arange(44, (180 / math.pi) * (math.atan(20 / 2)))
        eff_1 = 100 - (5 / 3) * tc_1
        eff_2 = 50 - (5 / 9) * tc_2

        ax.plot(tc_1, eff_1, color='green', linestyle='dashed', label='First Principles Line')
        ax.plot(tc_2, eff_2, color='green', linestyle='dashed')

    if air_gap == True:
        tc_1 = np.arange(0, 46)
        tc_2 = np.arange(44, (180 / math.pi) * (math.atan(20 / 2)))
        eff_1 = (5 / 9) * tc_1
        eff_2 = 50 - (5 / 9) * tc_2

        ax.plot(tc_1, eff_1, color='green', linestyle='dashed', label='First Principles Line')
        ax.plot(tc_2, eff_2, color='green', linestyle='dashed')


    ax.plot(tc, eff, color='blue', linestyle='solid', label='Numerical Simulation Line')
    ax.set(xlabel='Theta Critical (Degrees)', ylabel='Light Efficiency (eff)',
           title='Light Efficiency as a Function of Critical Angle')
    ax.set_ylim([0, 100])
    ax.set_xlim([0, 90])

    plt.text(70, 70, "l = 20.0")
    #plt.text(1.8, 85, "w = 3.0")
    plt.text(70, 65, "h = 2.0")
    plt.text(70, 60, "wp = " + str(wp))
    plt.text(70, 55, "hp = " + str(hp))
    plt.text(70, 50, "Detector = 1")
    plt.text(70, 45, "Airgap = " + str(air_gap))
    plt.text(40, 90, "Dotted line = First Principles")
    plt.text(40, 85, "Solid Line = Numerical Simulation")

    plt.show()
    plt.savefig("tc_vary.png")

def l_vary(wp, hp, air_gap):
    l = []
    eff = []

    for i in range(30):
        sim = Simulation(1.0 + i, 3.0, 2.0, 0.0, wp, hp, 1.6, 1.0, 1, air_gap)
        l.append(sim.l)
        eff.append(sim.run())

    fig, ax = plt.subplots()

    ax.plot(l, eff)

    ax.set(xlabel='Length of Bar (l)', ylabel='Light Efficiency (eff)',
               title='Light Efficiency as a Function of Length')
    ax.set_ylim([0, 100])
    ax.set_xlim([0, 30])

    #plt.text(22.5, 90, "w = 3.0")
    plt.text(22.5, 85, "h = 2.0")
    plt.text(22.5, 80, "wp = " + str(wp))
    plt.text(22.5, 75, "hp = " + str(hp))
    plt.text(22.5, 70, "n1 = 1.6")
    plt.text(22.5, 65, "Detector = 1")
    plt.text(22.5, 60, "Airgap = " + str(air_gap))

    plt.show()
    plt.savefig("l_vary.png")

def w_vary(wp, hp, air_gap):
    w = []
    eff = []

    for i in range(20):
        sim = Simulation(20.0, 2.0 + (i/10), 2.0, 0.0, wp, hp, 1.6, 1.0, 1, air_gap)
        w.append(sim.w)
        eff.append(sim.run())

    fig, ax = plt.subplots()

    ax.plot(w, eff)

    ax.set(xlabel='Width of Bar (w)', ylabel='Light Efficiency (eff)',
           title='Light Efficiency as a Function of Width')
    ax.set_ylim([0, 100])
    ax.set_xlim([2.0, 3.5])

    plt.text(3.1, 90, "l = 20.0")
    plt.text(3.1, 85, "h = 2.0")
    plt.text(3.1, 80, "wp = " + str(wp))
    plt.text(3.1, 75, "hp = " + str(hp))
    plt.text(3.1, 70, "n1 = 1.6")
    plt.text(3.1, 65, "Detector = 1")
    plt.text(3.1, 60, "Airgap = " + str(air_gap))

    plt.show()
    plt.savefig("w_vary.png")


def h_vary(wp, hp, air_gap):
    h = []
    eff = []

    for i in range(16):
        sim = Simulation(20.0, 3.0, 2.0 + (i/10), 0.0, wp, hp, 1.6, 1.0, 1, air_gap)
        h.append(sim.h)
        eff.append(sim.run())

    fig, ax = plt.subplots()

    ax.plot(h, eff)

    ax.set(xlabel='Height of Bar (h)', ylabel='Light Efficiency (eff)',
           title='Light Efficiency as a Function of Height')
    ax.set_ylim([0, 100])
    ax.set_xlim([2.0, 3.5])

    plt.text(3.1, 90, "l = 20.0")
    plt.text(3.1, 85, "w = 3.0")
    plt.text(3.1, 80, "wp = " + str(wp))
    plt.text(3.1, 75, "hp = " + str(hp))
    plt.text(3.1, 70, "n1 = 1.6")
    plt.text(3.1, 65, "Detector = 1")
    plt.text(3.1, 60, "Airgap = " + str(air_gap))

    plt.show()
    plt.savefig("h_vary.png")


def calculated_plot(airgap):
    if airgap == False:
        tc_1 = np.arange(0, 46)
        tc_2 = np.arange(44, (180/math.pi)*(math.atan(20/2)))
        eff_1 = 100 - (5 / 3) * tc_1
        eff_2 = 50 - (5 / 9) * tc_2
        eff_3 = 50 - (5 / 9) * (180 / math.pi) * (math.atan(20 / 2))

        fig, ax = plt.subplots()
        ax.set(xlabel='Theta Critical (Degrees)', ylabel='Light Efficiency (eff)',
               title='Light Efficiency as a Function of Critical Angle (Calculated)')

        ax.set_ylim([0, 100])
        ax.set_xlim([0, 90])
        ax.plot(tc_1, eff_1)
        ax.plot(tc_2, eff_2)
        ax.hlines(y=(50 - (5 / 9) * (180 / math.pi) * (math.atan(20 / 2))), xmin=(180 / math.pi) * (math.atan(20 / 2)),
                  xmax=90, color='r', linestyle='-')
        plt.show()
        plt.show()

    if airgap == True:
        tc_1 = np.arange(0, 46)
        tc_2 = np.arange(44, (180/math.pi)*(math.atan(20/2)))
        #tc_3 = np.arange((180/math.pi)*(math.atan(20/2)), 90)
        eff_1 = (5/9) * tc_1
        eff_2 = 50 - (5/9) * tc_2
        eff_3 = 50 - (5/9) * (180/math.pi)*(math.atan(20/2))

        fig, ax = plt.subplots()
        ax.set(xlabel='Theta Critical (Degrees)', ylabel='Light Efficiency (eff)',
               title='Light Efficiency as a Function of Critical Angle (Calculated)')

        ax.set_ylim([0, 100])
        ax.set_xlim([0, 90])
        ax.plot(tc_1, eff_1)
        ax.plot(tc_2, eff_2)
        #ax.plot(tc_3, eff_3)
        ax.hlines(y=(50 - (5/9) * (180/math.pi)*(math.atan(20/2))), xmin=(180/math.pi)*(math.atan(20/2)), xmax=90, color='r', linestyle='-')
        plt.show()


def photon_histogram(tc, air_gap):
    sim = Simulation(20.0, 3.0, 2.0, 0.0, 3.0, 2.0, 1 / math.sin(tc * (math.pi / 180)), 1.0, 1, air_gap)
    photon_data = sim.check()
    photon_pass = photon_data[0]
    photon_hit = photon_data[1]
    phi_select = photon_data[2]

    if 0 <= tc < 45:
        bins = [0, tc, 90-tc, 90, 90+tc, 180-tc, 180]
        plt.hist(photon_pass, bins = bins, edgecolor="black", label='Pass')
        plt.hist(photon_hit, bins=bins, edgecolor='black', label='Hit')
        plt.title("Theta Critical=" + str(tc) + "; Airgap is " + str(air_gap) + "; Phi=" + str(phi_select)
                  + "; bins=" + str(bins))
        plt.xlabel("Angles")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()
    if 45 <= tc <= 90:
        bins = [0, tc, 90, 180-tc, 180] #[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
        plt.hist(photon_pass, bins=bins, edgecolor="black", label='Pass')
        plt.hist(photon_hit, bins=bins, edgecolor='black', label='Hit')
        plt.title("Theta Critical=" + str(tc) + "; Airgap is " + str(air_gap) + "; Phi=" + str(phi_select)
                  + "; bins=" + str(bins))
        plt.xlabel("Angles")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()

def photon_pathlength():

    path_length = []
    eff = []

    for i in range(100):
        sim = Simulation(30.0, 2.0, 3.0, 0.0, 2.0, 3.0, 1.6, 1.0, 1, False, random.uniform(-math.pi/2, math.pi/2),
                         random.uniform(0, math.pi), random.uniform(-0.5, 0.5), random.uniform(-0.5, 0.5))
        path_length.append(sim.path_length())
        eff.append(sim.eff())

    bins = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    plt.hist(path_length, bins=bins, edgecolor="black", label='Pass')
    plt.title("Path Length Distribution")
    plt.xlabel("Lengths")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()
    #print(path_length)



#nominal dimensions of bar: l=20.0, w=3.0, h=2.0
#n1_vary(2.0, 1.0, False)
#tc_vary(3.0, 2.0, True)
#l_vary(3.0, 2.0, False)
#w_vary(2.0, 1.0, True)
#h_vary(2.0, 1.0, False)
#calculated_plot(False)
#photon_histogram(23, False)
photon_pathlength()
