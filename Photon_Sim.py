import math
import numpy as np
import random

class Simulation:

    def __init__(self, l, w, h, lp, wp, hp, n1, n2, detector, air_gap, phi_line,
                 theta_line, Xoy, Xoz):
        self.l = l # scintillator length: x
        self.w = w # scintillator width: y
        self.h = h # scintillator height: z
        self.lp = lp # ? length
        self.wp = wp # ? width
        self.hp = hp # ? height
        self.n1 = n1 # scintillator index of refraction
        self.n2 = n2 # air index of refraction
        self.detector = detector # integer 0-5 indicating which face if scintillator attaches to SiPM
        self.air_gap = air_gap # True if there is air between
        self.phi_line = phi_line
        self.theta_line = theta_line
        self.Xoy = Xoy # displacement from center of scintillator?
        self.Xoz = Xoz # displacement from center of scintillator?
        self.theta_critical = (math.asin(n2 / n1)) # minimum angle for TIR


    def random_line(self):

        # this seems to simulate random starting lines for a photon to initialize on rather than a path of an electron that this photons start at

        Tmin_x = -self.l / 2
        Tmax_x = self.l / 2
        Tmin_y = ((-self.w / 2) - self.Xoy) / math.tan(self.phi_line)
        Tmax_y = ((self.w / 2) - self.Xoy) / math.tan(self.phi_line)
        Tmin_z = ((-self.h / 2) - self.Xoz) * math.tan(self.theta_line)
        Tmax_z = ((self.h / 2) - self.Xoz) * math.tan(self.theta_line)
        # checks the T boundaries of each x, y, and z based on the dimensions of the box, l, w, and h

        self.Ts = [Tmin_x, Tmax_x, Tmin_y, Tmax_y, Tmin_z, Tmax_z]
        self.Ts.sort()
        # sorts out all possible T values; the two middle values after sorting are the Tmin and Tmax that would satisfy
        # the range of T's of the generated random line

        for i in range(2, 4):
            if ((-self.l / 2 <= self.Ts[i] <= self.l / 2) and (
                    -self.w / 2 <= self.Ts[i] * math.tan(self.phi_line) + self.Xoy <= self.w / 2) and
                    (-self.h / 2 <= (self.Ts[i]/math.tan(self.theta_line)) + self.Xoz <= self.h / 2)):
                self.length = (self.Ts[3] - self.Ts[2]) * math.sqrt(
                    1 + (math.tan(self.phi_line)) ** 2 + 1/(math.tan(self.theta_line)) ** 2)
            else:
                #print("T is not defined")
                return False


    def photon(self, V, Ro, rec=0):

        if rec > 900:
            return False

        for i in range(1):
            if V[i] > 0.:
                x = self.l/2
                t = (x - Ro[i]) / V[i]
                y = Ro[i+1] + V[i+1] * t
                z = Ro[i+2] + V[i+2] * t
                R = np.array([x, y, z])  #finds the coordinates where the photon hits the x = l/2 plane
            if V[i] < 0.:
                x = -self.l/2
                t = (x - Ro[i]) / V[i]
                y = Ro[i+1] + V[i+1] * t
                z = Ro[i+2] + V[i+2] * t
                R = np.array([x, y, z])  #finds the coordinates where the photon hits the x = -l/2 plane
            if (-self.l/2 <= R[i] <= self.l/2) and (-self.w/2 <= R[i+1] <= self.w/2) and (-self.h/2 <= R[i+2] <= self.h/2): #checks to see if any of those two points are within the boundaries of the box
                theta_i = (math.acos(abs(V[i]) / math.sqrt(V[i] ** 2 + V[i+1] ** 2 + V[i+2] ** 2)))
                if theta_i > self.theta_critical:
                    if self.air_gap:
                        V[i] *= -1
                        Ro = R
                        return self.photon(V, Ro, rec+1)
                    if not(self.air_gap):
                        if (-self.wp <= R[i+1] <= self.wp and -self.hp <= R[i+2] <= self.hp) and \
                        ((self.detector == 1 and x == self.l / 2) or (self.detector == 2 and x == -self.l / 2)):
                            return True
                        else:
                            V[i] *= -1
                            Ro = R
                            return self.photon(V, Ro, rec+1)
                else:
                    theta_t = math.asin(self.n1 * math.sin(theta_i))
                    r_perp = (self.n1 * math.cos(theta_i) - self.n2 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n2 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n2 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 101)
                    Min = min(Transmittance, Reflectance)
                    if self.air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro, rec+1)
                            if Min == Transmittance:
                                if (-self.wp <= R[i + 1] <= self.wp and -self.hp <= R[i + 2] <= self.hp) and \
                                ((self.detector == 1 and x == self.l / 2) or (self.detector == 2 and x == -self.l / 2)):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (-self.wp <= R[i + 1] <= self.wp and -self.hp <= R[i + 2] <= self.hp) and \
                                ((self.detector == 1 and x == self.l / 2) or (self.detector == 2 and x == -self.l / 2)):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro, rec+1)
                    if not self.air_gap:
                        if (-self.wp <= R[i + 1] <= self.wp and -self.hp <= R[i + 2] <= self.hp) and \
                        ((self.detector == 1 and x == self.l / 2) or (self.detector == 2 and x == -self.l / 2)):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro, rec+1)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro, rec+1)


        for i in range(1,2):
            if V[i] > 0.:
                y = self.w/2
                t = (y - Ro[i]) / V[i]
                x = Ro[i-1] + V[i-1] * t
                z = Ro[i+1] + V[i+1] * t
                R = np.array([x, y, z])
            if V[i] < 0.:
                y = -self.w/2
                t = (y - Ro[i]) / V[i]
                x = Ro[i-1] + V[i-1] * t
                z = Ro[i+1] + V[i+1] * t
                R = np.array([x, y, z])
            if (-self.l/2 <= R[i-1] <= self.l/2) and (-self.w/2 <= R[i] <= self.w/2) and (-self.h/2 <= R[i+1] <= self.h/2):
                theta_i = (math.acos(abs(V[i]) / math.sqrt(V[i-1] ** 2 + V[i] ** 2 + V[i+1] ** 2)))
                if theta_i > self.theta_critical:
                    if self.air_gap:
                        V[i] *= -1
                        Ro = R
                        return self.photon(V, Ro, rec+1)
                    if not(self.air_gap):
                        if (self.detector == 3 and y == self.w / 2) or (self.detector == 4 and y == -self.w / 2):
                            return True
                        else:
                            V[i] *= -1
                            Ro = R
                            return self.photon(V, Ro, rec+1)
                else:
                    theta_t = math.asin(self.n1 * math.sin(theta_i))
                    r_perp = (self.n1 * math.cos(theta_i) - self.n2 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n2 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n2 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 100)
                    Min = min(Transmittance, Reflectance)
                    if self.air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro, rec+1)
                            if Min == Transmittance:
                                if (self.detector == 3 and y == self.w / 2) or (self.detector == 4 and y == -self.w / 2):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (self.detector == 3 and y == self.w / 2 ) or (self.detector == 4 and y == -self.w / 2):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro, rec+1)
                    if not self.air_gap:
                        if (self.detector == 3 and y == self.w / 2) or (self.detector == 4 and y == -self.w / 2):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro, rec+1)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro, rec+1)

        for i in range(2,3):
            if V[i] > 0.:
                z = self.h/2
                t = (z - Ro[i]) / V[i]
                x = Ro[i-2] + V[i-2] * t
                y = Ro[i-1] + V[i-1] * t
                R = np.array([x, y, z])
            if V[i] < 0.:
                z = -self.h/2
                t = (z - Ro[i]) / V[i]
                x = Ro[i-2] + V[i-2] * t
                y = Ro[i-1] + V[i-1] * t
                R = np.array([x, y, z])
            if (-self.l/2 <= R[i-2] <= self.l/2) and (-self.w/2 <= R[i-1] <= self.w/2) and (-self.h/2 <= R[i] <= self.h/2):
                theta_i = (math.acos(abs(V[i]) / math.sqrt(V[i-2] ** 2 + V[i-1] ** 2 + V[i] ** 2)))
                if theta_i > self.theta_critical:
                    if self.air_gap:
                        V[i] *= -1
                        Ro = R
                        return self.photon(V, Ro, rec+1)
                    if not(self.air_gap):
                        if (self.detector == 5 and z == self.h / 2) or (self.detector == 6 and z == -self.h / 2):
                            return True
                        else:
                            V[i] *= -1
                            Ro = R
                            return self.photon(V, Ro, rec+1)
                else:
                    theta_t = math.asin(self.n1 * math.sin(theta_i))
                    r_perp = (self.n1 * math.cos(theta_i) - self.n2 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n2 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n2 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) * 100
                    Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 101)
                    Min = min(Transmittance, Reflectance)
                    if self.air_gap:
                        if 0 <= select_path <= Min:
                            if Min == Reflectance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro, rec+1)
                            if Min == Transmittance:
                                if (self.detector == 5 and z == self.h / 2) or (self.detector == 6 and z == -self.h / 2):
                                    return True
                                else:
                                    return False
                        else:
                            if Min == Reflectance:
                                if (self.detector == 5 and z == self.h / 2) or (self.detector == 6 and z == -self.h / 2):
                                    return True
                                else:
                                    return False
                            if Min == Transmittance:
                                V[i] *= -1
                                Ro = R
                                return self.photon(V, Ro, rec+1)
                    if not self.air_gap:
                        if (self.detector == 5 and z == self.h / 2) or (self.detector == 6 and z == -self.h / 2):
                            return True
                        else:
                            if 0 <= select_path <= Min:
                                if Min == Reflectance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro, rec+1)
                                if Min == Transmittance:
                                    return False
                            else:
                                if Min == Reflectance:
                                    return False
                                if Min == Transmittance:
                                    V[i] *= -1
                                    Ro = R
                                    return self.photon(V, Ro, rec+1)

    def run(self, detected_photon=0):
        total = 1000
        self.random_line()

        self.photon_pass_phi_0 = []
        self.photon_hit_phi_0 = []
        self.photon_pass_phi_pi = []
        self.photon_hit_phi_pi = []

        for n in range(total):
            if self.random_line() == False:
                break
            phi_initial =  random.choice([0, math.pi]) #random.uniform(0, 2 * math.pi)
            theta_initial = random.uniform(0, math.pi)  # math.acos(random.uniform(1., -1))   #(90) * math.pi/180
            Vx = math.sin(theta_initial) * math.cos(phi_initial)
            Vy = math.sin(theta_initial) * math.sin(phi_initial)
            Vz = math.cos(theta_initial)
            V = np.array([Vx, Vy, Vz])
            # generates random direction (V) based on random thetas of phi and theta

            self.T = random.uniform(self.Ts[2], self.Ts[3])
            Rox = self.T  # x(t)
            Roy = (self.T - self.Xoy) * math.tan(self.phi_line)  # y(t)
            Roz = (self.T - self.Xoz) * math.tan((math.pi / 2) - self.theta_line)  # z(t)
            Ro = np.array([Rox, Roy, Roz])
            # generates random point (Ro) on a line based on the angle phi and theta of that line

            if phi_initial == 0:
                self.photon_pass_phi_0.append((180 / math.pi) * theta_initial)
                if self.photon(V, Ro) == True:
                    self.photon_hit_phi_0.append((180 / math.pi) * theta_initial)
            if phi_initial == math.pi:
                self.photon_pass_phi_pi.append((180 / math.pi) * theta_initial)
                if self.photon(V, Ro) == True:
                    self.photon_hit_phi_pi.append((180 / math.pi) * theta_initial)

            if self.photon(V, Ro) == True:
                   detected_photon += 1

        self.efficiency = (detected_photon * 100) / total


    def check(self):
        self.run()
        choose_phi = int(input("What phi angle?"))
        if choose_phi == 0:
            return self.photon_pass_phi_0, self.photon_hit_phi_0, choose_phi
        if choose_phi == 180:
            return self.photon_pass_phi_pi, self.photon_hit_phi_pi, choose_phi

    def path_length(self):
        self.run()
        if self.random_line() == False:
            return 0
        else:
            #print(self.length)
            return self.length

    def eff(self):
        self.run()
        return self.efficiency



#sim = Simulation(l, w, h, lp, wp, hp, n1, n2, detector, air_gap)
sim = Simulation(2.0, 30.0, 3.0, 0.0, 2.0, 3.0, 1.58, 1.0, 1, False, math.pi/4, math.pi/4, 0, 0)

sim.run()

print(f'Efficiency: {sim.eff()}%')
print(f'Path length: {sim.path_length()}')