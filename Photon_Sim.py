import math
import numpy as np
import random

class Simulation:

    def __init__(self, l, w, h, lp, wp, hp, n1, n2, n3, phi_line=math.pi/4,
                 theta_line=math.pi/4, detector=3, air_gap=False, Xoy=0, Xoz=0, iterations=1000, nplastic=1.502):
        self.l = l # scintillator length: x 2?
        self.w = w # scintillator width: y 30?
        self.h = h # scintillator height: z 3?
        self.lp = lp # SiPM window length
        self.wp = wp # SiPM window width
        self.hp = hp # SiPM window height
        self.n1 = n1 # scintillator index of refraction
        self.n2 = n2 # air index of refraction
        self.n3 = n3 # SiPM index of refraction
        self.detector = [False, False, False, False, False, False] # iterable boolean list indicating which face of the
        # scintillator attaches to the SiPM window
        self.detector[detector] = True
        self.air_gap = air_gap # True if there is air between
        self.phi_line = phi_line # phi of old random_line method
        self.theta_line = theta_line # theta of old random_line method
        self.Xoy = Xoy # displacement from center of scintillator?
        self.Xoz = Xoz # displacement from center of scintillator?
        self.theta_critical = (math.asin(n2 / n1)) # minimum angle for TIR
        self.iterations = iterations
        self.nback = nplastic
        #self.theta_pass = (math.asin(nplastic / n1))
        self.theta_detect = (math.asin(n3 / n1))
        self.back = False # allows or disallows backflow
        self.theta_back = math.pi / 2
        self.lwhb = [0, 0, 0] # backflow window length, width, and height

    def ray_trace(self, V, Ro, rec=0, length=0):

        # This method takes simulation inputs and photon position and velocity NumPy Arrays to recursively raytrace
        # which wall it will hit and at what angle until it is reasonably absorbed or escapes, resulting in a False,
        # or it is detected, resulting in a True.
        # These are returned as the first element in a list so that more information can be returned if the photon
        # must pass through multiple objects
        # Returns: [boolean: true if photon passes desired window,
        #           float: length photon traveled,
        #           array of floats(if passes): position of passage,
        #           array of floats(if passes): exit velocity vector,
        #           boolean: true if backflows to previous stage]
        # Detector = 2 (or 3) are correct for the current dimensionality inputs

        if length > 3800 or rec > 900: # !rewrite later to take attenuation length into account (attenuation length: 380 cm)
            #print(f'Absorbed in {rec}')
            #print(length)
            return [False, length, False]

        dims = [self.l, self.w, self.h]
        window = [self.lp, self.wp, self.hp]  # allows easier iterating across dimensions

        for i in range(3): # checks each wall of the scintillator until it finds the one that the photon will hit
            # i equalling 0 in this loop makes this section check the x component of V, and so on.
            if V[i] == 0:
                continue

            R = np.zeros(3)
            R[i] = math.copysign(dims[i]/2, V[i])
            t = (R[i] - Ro[i]) / V[i]
            R[(i+1) % 3] = Ro[(i+1) % 3] + V[(i+1) % 3] * t
            R[(i+2) % 3] = Ro[(i+2) % 3] + V[(i+2) % 3] * t
            #finds the coordinates where the photon hits the plane of each of the scintillators walls

            if (np.abs(R[i]) <= dims[i]/2) and (np.abs(R[(i+1) % 3]) <= dims[(i+1) % 3]/2) and (np.abs(R[(i+2) % 3]) <= dims[(i+2) % 3]/2): # checks to see if the new point is within the boundaries of the box
                theta_i = (math.acos(abs(V[i]) / math.sqrt(V[(i+1) % 3] ** 2 + V[i] ** 2 + V[(i+2) % 3] ** 2))) # extracts angle the photon intersects with the wall
                #print(theta_i)
                length += np.linalg.norm(R - Ro)
                #print(length)

                if np.abs(R[(i + 1) % 3]) <= window[(i + 1) % 3] / 2 and np.abs(R[(i + 2) % 3]) <= window[(i + 2) % 3] / 2 \
                    and ((self.detector[i * 2] and R[i] == dims[i] / 2) or (self.detector[(i * 2) + 1] and R[i] == -dims[i] / 2)):  # checks that the photon could hit the detector at this intersection point
                    if theta_i > self.theta_detect:
                        # Immediately returns false because a TIR bounce on the necessary passage geometrically
                        # disallows a photon from ever crossing this threshold for rectangular geometry
                        return [False, length, False]
                    theta_t = math.asin((self.n1 / self.n3) * math.sin(theta_i))  # transmission angle
                    r_perp = (self.n1 * math.cos(theta_i) - self.n3 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n3 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n3 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n3 * math.cos(theta_i))  # Fresnel's equations
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2)  # calculates average reflectance regardless of polarization and converts into decimal chance of reflecting
                    select_path = random.uniform(0, 1)
                    if select_path <= Reflectance:
                        V[i] *= -1
                        # print('rare bounce')
                        return self.ray_trace(V, R, rec + 1, length)
                    # calculate new V
                    phi = np.arctan(V[(i + 2) % 3] / V[(i + 1) % 3])
                    V[(i + 1) % 3] = math.copysign(np.sin(theta_t) * np.cos(phi), V[(i + 1) % 3])
                    V[(i + 2) % 3] = math.copysign(np.sin(theta_t) * np.sin(phi), V[(i + 2) % 3])
                    V[i] = math.copysign(np.cos(theta_t), V[i])
                    return [True, length, R, V, False]

                if np.abs(R[(i + 1) % 3]) <= self.lwhb[(i + 1) % 3] / 2 and np.abs(R[(i + 2) % 3]) <= self.lwhb[(i + 2) % 3] / 2 \
                    and self.back and ((self.detector[i * 2] and R[i] == -dims[i] / 2) or (self.detector[(i * 2) + 1] and R[i] == dims[i] / 2)):  # checks if the photon could backflow at this intersection point
                    #print(f'back: {self.lwhb[(i + 1) % 3] / 2} by {self.lwhb[(i + 2) % 3] / 2} \n ')
                    if theta_i > self.theta_back:
                        # print('TIR bounce')
                        V[i] *= -1
                        return self.ray_trace(V, R, rec + 1, length)
                    theta_t = math.asin((self.n1 / self.nback) * math.sin(theta_i))  # transmission angle
                    r_perp = (self.n1 * math.cos(theta_i) - self.nback * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.nback * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.nback * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.nback * math.cos(theta_i))  # Fresnel's equations
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2)  # calculates average reflectance regardless of polarization and converts into decimal chance of reflecting
                    select_path = random.uniform(0, 1)
                    if select_path <= Reflectance:
                        V[i] *= -1
                        # print('rare bounce')
                        return self.ray_trace(V, R, rec + 1, length)
                    # calculate new V
                    phi = np.arctan(V[(i + 2) % 3] / V[(i + 1) % 3])
                    V[(i + 1) % 3] = math.copysign(np.sin(theta_t) * np.cos(phi), V[(i + 1) % 3])
                    V[(i + 2) % 3] = math.copysign(np.sin(theta_t) * np.sin(phi), V[(i + 2) % 3])
                    V[i] = math.copysign(np.cos(theta_t), V[i])
                    #print(V)
                    return [False, length, R, V, True]

                if theta_i > self.theta_critical: # avoids extra computation for case of TIR

                    # print('TIR bounce')
                    V[i] *= -1
                    return self.ray_trace(V, R, rec+1, length)
                else:
                    theta_t = math.asin((self.n1 / self.n2) * math.sin(theta_i)) # transmission angle
                    r_perp = (self.n1 * math.cos(theta_i) - self.n2 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n2 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n2 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n2 * math.cos(theta_i)) # Fresnel's equations
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) # calculates average reflectance regardless of polarization and converts into decimal chance of reflecting
                    #Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 1)

                    if select_path <= Reflectance:
                        V[i] *= -1
                        # print('bounce')
                        return self.ray_trace(V, R, rec + 1, length)
                    # print('escape')
                    return [False, length, False]

        raise Exception(f'Photon tunneled out of sim, look for bugs \n R: {Ro} \n V: {V}, \n Dims: {dims}')

    def random_three_vector(self):
        # written by Majd Ghrear in a previous project to isotropically generate direction vectors
        phi = np.random.uniform() * 2 * np.pi

        costheta = 2.0 * np.random.uniform() - 1.0
        theta = np.arccos(costheta)

        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)

        return np.array([x, y, z]), theta, phi

    def random_test(self, plastic=False):

        # generates self.iterations number of photons with random position and velocities withing the scintillator and returns the fraction that are detected

        count = 0
        dims = np.array([self.l, self.w, self.h])
        hits = []
        misses = []
        for i in range(self.iterations):
            Ro = np.random.uniform(low=-1.0, high=1.0, size=3) * dims / 2
            #Ro = self.random_three_vector()[0] * dims / 2
            Vo = self.random_three_vector()[0]
            detection = self.ray_trace(Vo, Ro)
            if detection[0]:
                if plastic:
                    Ro = detection[1]
                else:
                    hits.append(detection[1:2])
                    count += 1
            else:
                misses.append(detection[1:2])
        return count / self.iterations, np.array(hits), np.array(misses)

    def input_test(self, y, z, n=None):

        # generates self.iterations number of photons with random position and velocities withing the scintillator and returns the fraction that are detected

        if n is None:
            n = self.iterations
        count = 0
        dims = np.array([self.l, self.w, self.h])
        for i in range(n):
            Ro = np.array([np.random.uniform(low=-1.0, high=1.0) * dims[0] / 2, y, z])
            Vo = self.random_three_vector()[0]
            detection = self.ray_trace(Vo, Ro)
            if detection[0]:
                count += 1
        return count / n

    def run(self, y, z, dimensions, *args, n=None):

        # Generates self.iterations number of photons with random positions on an input line and random velocities
        # within the scintillator.
        # Takes coordinates, dimensions for the intermediate portions, and arguments for the intermediate indices of
        # refraction and returns the fraction that are detected.
        # Assumes detector == 2 or 3

        if n is None:
            n = self.iterations
        count = 0
        dims = np.concatenate((np.array([[self.l, self.w, self.h]]), dimensions, np.array([[self.lp, self.wp, self.hp]])))
        r_indices = [self.n1] + list(args) + [self.n3] # list indices of refraction to allow iteration through the
        # mediums a photon must transit and allows restoration after run
        #print(dims)
        for i in range(n):
            Ro = np.array([np.random.uniform(low=-1.0, high=1.0) * dims[0, 0] / 2, y, z])
            Vo = self.random_three_vector()[0]
            j = 0 # tracks which stage the photon is in
            length = 0
            self.theta_back = math.pi / 2
            while j < len(r_indices) - 1:
                #print(j)
                # iterating through index of refraction
                self.n1 = r_indices[j]
                self.n3 = r_indices[j + 1]
                self.theta_critical = (math.asin(self.n2 / self.n1))
                try:
                    self.theta_detect = (math.asin(self.n3 / self.n1))
                except ValueError:
                    self.theta_detect = math.pi / 2

                if j > 0: # sets backflow values when backflow is possible and turns backflow off when it is not
                    self.back = True
                    self.nback = r_indices[j - 1]
                    try:
                        self.theta_back = (math.asin(self.nback / self.n1))
                    except ValueError:
                        self.theta_back = math.pi / 2
                    self.lwhb = [dims[j - 1, 0], dims[j - 1, 1], dims[j - 1, 2]]
                else:
                    self.back = False

                # iterating through dimensions
                self.l, self.w, self.h = dims[j, 0], dims[j, 1],dims[j, 2]
                self.lp, self.wp, self.hp = dims[j + 1, 0], dims[j + 1, 1], dims[j + 1, 2]

                detection = self.ray_trace(Vo, Ro, length=length)
                if detection[0] or detection[-1]:
                    length += detection[1]
                    if detection[-1]: # backflow condition
                        j -= 2
                        Ro, Vo = detection[2], detection[3]
                        Ro[1] = dims[j + 1, 1] / 2
                        # print('backflow')
                        # print(Ro)
                        # print(Vo)
                    elif j == len(r_indices) - 2: # detection condition
                        count += 1
                    else: # passage to next stage condition
                        #print(detection[2:4])
                        Ro, Vo = detection[2], detection[3]
                        Ro[1] = - dims[j + 1, 1] / 2
                        #print(Ro)
                        #print(Vo)
                else:
                    break
                j += 1

        # reset all instance values so run method can be reused
        self.l, self.w, self.h = dims[0, 0], dims[0, 1], dims[0, 2]
        self.lp, self.wp, self.hp = dims[-1, 0], dims[-1, 1], dims[-1, 2]
        self.n1 = r_indices[0]
        self.n3 = r_indices[-1]
        self.theta_critical = (math.asin(self.n2 / self.n1))
        self.back = False
        return count / n




#sim = Simulation(l, w, h, lp, wp, hp, n1, n2, phi_line, theta_line)
sim = Simulation(2.0, 30.0, 3.0, 2.0, 30.0, 2.0, 1.58, 1.0, 1.55, detector=2)

#sim.run()
#print(f'Efficiency: {sim.efficiency}%')
#sim.new_line()
#print(f'Path length: {sim.length}')
#print(f'Path length new: {sim.path_length()}') # currently unrelated to previous run

# V = np.array([0, 1, 2])
# Ro = np.array([0, 0, 0])
# print(sim.theta_critical)
# if sim.ray_trace(V, Ro):
#     print('Detected')
# else:
#     print('Lost')

#print(f'Detected {sim.random_test()[0] * 100}%')
#print(f'Detected {sim.input_test(0, 0) * 100}%')


dimensions = np.array([[2.0, 0.125, 3.0], [2.0, 54.86, 3.0], [100.0, 0.1, 100.0]])
print(sim.run(0, 0, dimensions, 1.57, 1.502, 1.0))
