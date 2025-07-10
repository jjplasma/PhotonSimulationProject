import math
import numpy as np
import random

class Simulation:

    def __init__(self, l, w, h, lp, wp, hp, n1, n2, n3, phi_line,
                 theta_line, detector=3, air_gap=False, Xoy=0, Xoz=0, iterations=1000, nplastic=1.502):
        self.l = l # scintillator length: x 2?
        self.w = w # scintillator width: y 30?
        self.h = h # scintillator height: z 3?
        self.lp = lp # SiPM window length
        self.wp = wp # SiPM window width
        self.hp = hp # SiPM window height
        self.n1 = n1 # scintillator index of refraction
        self.n2 = n2 # air index of refraction
        self.n3 = n3 # SiPM index of refraction
        self.detector = [False, False, False, False, False, False] # iterable boolean list indicating which face if scintillator attaches to SiPM
        self.detector[detector] = True
        self.air_gap = air_gap # True if there is air between
        self.phi_line = phi_line
        self.theta_line = theta_line
        self.Xoy = Xoy # displacement from center of scintillator?
        self.Xoz = Xoz # displacement from center of scintillator?
        self.theta_critical = (math.asin(n2 / n1)) # minimum angle for TIR
        self.iterations = iterations
        self.nplastic = nplastic
        self.theta_pass = (math.asin(nplastic / n1))
        #self.theta_detect = (math.asin(n3 / n1))


    def random_line(self):

        # this seems to simulate random starting lines for a photon to initialize on rather than a path of an electron that the photons start at
        # no randomness occurs within this method, rather it updates self.length to the most recent attribute values or returns false if they lead to an impossible path

        Tmin_x = -self.l / 2
        Tmax_x = self.l / 2 # x direction originally was longest direction, but was switched to beam direction since it works differently here and that's the only direction that should work differently
        #print(self.phi_line)
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

    def new_line(self, phi=None, theta=None):

        # Adjustments to random_line

        # Define boundaries of the scintillator based on its width and height
        Tmin_x = -self.l / 2
        Tmax_x = self.l / 2

        # Initialize the boundaries for Y and Z
        ymin = -self.w / 2
        ymax = self.w / 2
        zmin = -self.h / 2
        zmax = self.h / 2

        # Based on the angle of theta, determine the path potential in Z:
        if self.theta_line == 0:  # Path is straight up in Z
            Tmin_z = zmin
            Tmax_z = zmax
            Tmin_y = ymin  # Y does not travel; stays fixed
            Tmax_y = ymax

        elif self.theta_line == math.pi:  # Path is straight down in Z
            Tmin_z = zmax  # We're taking the top boundary of Z
            Tmax_z = zmin  # We're going down to the bottom boundary
            Tmin_y = ymin
            Tmax_y = ymax

        else:  # For angles not strictly vertical
            # Calculate Tmin and Tmax for Y based on the phi angle
            if self.phi_line != 0:  # Avoid tangent if phi is 0
                Tmin_y = ((-self.w / 2) - self.Xoy) / math.tan(self.phi_line)
                Tmax_y = ((self.w / 2) - self.Xoy) / math.tan(self.phi_line)
            else:  # If phi is directly 0, Y component does not change
                Tmin_y = ymin
                Tmax_y = ymax

            # Z calculations now take into account the relationship to theta
            Tmin_z = ((-self.h / 2) - self.Xoz) * math.tan(self.theta_line)
            Tmax_z = ((self.h / 2) - self.Xoz) * math.tan(self.theta_line)

        # Store all T values and sort them
        self.Ts = [Tmin_x, Tmax_x, Tmin_y, Tmax_y, Tmin_z, Tmax_z]
        self.Ts.sort()  # Ensure T values are ordered

        # Initialize a validity flag
        valid = False

        # Check for valid indices for Y and Z
        for i in range(len(self.Ts)):
            # Check boundaries against scintillator dimensions appropriately
            if (-self.l / 2 <= self.Ts[i] <= self.l / 2):
                if i == 2 or i == 3:  # Only check Y boundaries
                    # Calculate path length based on Y values
                    self.length = min(Tmax_y, ymax) - max(Tmin_y, ymin)
                    valid = True
                    break
                elif i == 4 or i == 5:  # Only check Z boundaries
                    # Calculate Z length
                    self.length = min(Tmax_z, zmax) - max(Tmin_z, zmin)
                    valid = True
                    break

        return valid if valid else False


    def photon(self, V, Ro, rec=0):

        # This method takes simulation inputs and photon position and velocity to recursively raytrace
        # which wall it will hit and at what angle until it is reasonably absorbed or escapes, resulting in a False,
        # or it is detected, resulting in a True
        # It seems detector = 3 or 4 are correct for the current dimensionality inputs

        if rec > 900:
            #print('Absorbed')
            return False

        for i in range(1): # why are we looping once? seems to just set i = 0 but this could be done less confusingly
            if V[i] > 0.: # i being 0 in this loop makes this section check the x component of V
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
            if (-self.l/2 <= R[i] <= self.l/2) and (-self.w/2 <= R[i+1] <= self.w/2) and (-self.h/2 <= R[i+2] <= self.h/2): # checks to see if any of those two points are within the boundaries of the box
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


        for i in range(1,2): # i being 1 in this loop makes this section check the y component of V
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

        for i in range(2,3): # i being 2 in this loop makes this section check the z component of V
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

    def ray_trace(self, V, Ro, rec=0):

        # This method takes simulation inputs and photon position and velocity NumPy Arrays to recursively raytrace
        # which wall it will hit and at what angle until it is reasonably absorbed or escapes, resulting in a False,
        # or it is detected, resulting in a True
        # Detector = 2 (or 3) are correct for the current dimensionality inputs

        if rec > 900: # !rewrite later to take attenuation length into account
            print('Absorbed')
            return False

        for i in range(3): # checks each wall of the scintillator until it finds the one that the photon will hit
            # i equalling 0 in this loop makes this section check the x component of V, and so on.
            if V[i] == 0:
                continue

            dims = [self.l, self.w, self.h]
            window = [self.lp, self.wp, self.hp]
            R = np.zeros(3)
            R[i] = math.copysign(dims[i]/2, V[i])
            t = (R[i] - Ro[i]) / V[i]
            R[(i+1) % 3] = Ro[(i+1) % 3] + V[(i+1) % 3] * t
            R[(i+2) % 3] = Ro[(i+2) % 3] + V[(i+2) % 3] * t
            #finds the coordinates where the photon hits the plane of each of the scintillators walls

            if (np.abs(R[i]) <= dims[i]/2) and (np.abs(R[(i+1) % 3]) <= dims[(i+1) % 3]/2) and (np.abs(R[(i+2) % 3]) <= dims[(i+2) % 3]/2): # checks to see if any of those two points are within the boundaries of the box
                theta_i = (math.acos(abs(V[i]) / math.sqrt(V[(i+1) % 3] ** 2 + V[i] ** 2 + V[(i+2) % 3] ** 2))) # extracts angle the photon intersects with the wall
                # print(theta_i)
                if theta_i > self.theta_critical: # avoids extra computation for case of TIR
                    # if self.air_gap: # !this airgap detection will need adjusting for new plastic light pipe case
                    #     V[i] *= -1
                    #     return self.ray_trace(V, R, rec+1)
                    if np.abs(R[(i+1) % 3]) <= window[(i+1) % 3]/2 and np.abs(R[(i+1) % 3]) <= window[(i+2) % 3]/2 and \
                    ((self.detector[i * 2] and R[i] == dims[i] / 2) or (self.detector[(i * 2) + 1] and R[i] == -dims[i] / 2)): # checks that the photon could hit the detector at this intersection point
                        # Current assumption is that with the lower IOR difference between the SiPM and Scintillator, all photons will transmit
                        # print(f'if {theta_i} < 1.37')
                        # print('detected! \n')
                        return True
                    else:
                        # print('TIR bounce')
                        V[i] *= -1
                        return self.ray_trace(V, R, rec+1)
                else:
                    theta_t = math.asin(self.n1 * math.sin(theta_i)) # transmission angle
                    r_perp = (self.n1 * math.cos(theta_i) - self.n2 * math.cos(theta_t)) / (self.n1 * math.cos(theta_i) + self.n2 * math.cos(theta_t))
                    r_para = (self.n1 * math.cos(theta_t) - self.n2 * math.cos(theta_i)) / (self.n1 * math.cos(theta_t) + self.n2 * math.cos(theta_i))
                    r_ave = (abs(r_perp) + abs(r_para)) / 2
                    Reflectance = (r_ave ** 2) # calculates average reflectance regardless of polarization and converts into decimal chance of reflecting
                    #Transmittance = 100 - Reflectance
                    select_path = random.uniform(0, 1)


                    if np.abs(R[(i+1) % 3]) <= window[(i+1) % 3]/2 and np.abs(R[(i+1) % 3]) <= window[(i+2) % 3]/2 and \
                    ((self.detector[i * 2] and R[i] == dims[i] / 2) or (self.detector[(i * 2) + 1] and R[i] == -dims[i] / 2)): # checks that the photon could hit the detector at this intersection point
                        # print('detected! \n')
                        return True
                    if select_path <= Reflectance:
                        V[i] *= -1
                        # print('bounce')
                        return self.ray_trace(V, R, rec + 1)
                    # print('escape')
                    return False

        raise Exception('Photon tunneled out of sim, look for bugs')

    def random_test(self):

        # generates self.iterations number of photons with random position and velocities withing the scintillator and returns the fraction that are detected

        count = 0
        dims = np.array([self.l, self.w, self.h])
        for n in range(self.iterations):
            Ro = np.random.uniform(low=-1.0, high=1.0, size=3) * dims / 2
            Vo = np.random.uniform(low=-1.0, high=1.0, size=3)
            if self.ray_trace(Vo, Ro):
                count += 1
        return count / self.iterations


    def run(self, detected_photon=0):

        # runs the simulation of default 1000 photons emitting in the scintillator

        total = self.iterations
        # self.random_line() # what does this return? commented out because I believe this does nothing !update: this update self.length, but is redundant since it is called later

        self.photon_pass_phi_0 = []
        self.photon_hit_phi_0 = []
        self.photon_pass_phi_pi = []
        self.photon_hit_phi_pi = []#  Lists of angles photons travel at in certain conditions, unclear how this is used?

        for n in range(total):
            if self.random_line() == False: # maybe checks if the displacement inputs create a line outside the scintillator?
                break
            phi_initial =  random.choice([0, math.pi]) #random.uniform(0, 2 * math.pi)
            theta_initial = random.uniform(0, math.pi)  # math.acos(random.uniform(1., -1))   #(90) * math.pi/180
            Vx = math.sin(theta_initial) * math.cos(phi_initial)
            Vy = math.sin(theta_initial) * math.sin(phi_initial)
            Vz = math.cos(theta_initial)
            V = np.array([Vx, Vy, Vz])
            # generates random direction (V) based on random thetas of phi and theta

            # test code
            #self.phi_line = phi_initial
            #print(phi_initial)
            #print(self.phi_line)
            #self.theta_line = theta_initial
            #self.random_line()

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

        #print(self.photon_pass_phi_0)
        #print(self.photon_hit_phi_0)
        #print(self.photon_pass_phi_pi)
        #print(self.photon_hit_phi_pi)


    def check(self):
        self.run()
        choose_phi = int(input("What phi angle?"))
        if choose_phi == 0:
            return self.photon_pass_phi_0, self.photon_hit_phi_0, choose_phi
        if choose_phi == 180:
            return self.photon_pass_phi_pi, self.photon_hit_phi_pi, choose_phi

    def path_length(self):
        self.run()
        if not self.random_line():
            return 0
        else:
            #print(self.length)
            return self.length

    def eff(self):
        self.run()
        return self.efficiency



#sim = Simulation(l, w, h, lp, wp, hp, n1, n2, phi_line, theta_line)
sim = Simulation(2.0, 30.0, 3.0, 2.0, 30.0, 3.0, 1.58, 1.0, 1.55, math.pi/4, math.pi/2, detector=2)

#sim.run()
#print(f'Efficiency: {sim.efficiency}%')
#sim.new_line()
#print(f'Path length: {sim.length}')
#print(f'Path length new: {sim.path_length()}') # currently unrelated to previous run

# V = np.array([0, 1, 2]) #suspected bug with bounces on the top or bottom wall
# Ro = np.array([0, 0, 0])
# print(sim.theta_critical)
# if sim.ray_trace(V, Ro):
#     print('Detected')
# else:
#     print('Lost')

print(f'Detected {sim.random_test() * 100}%')