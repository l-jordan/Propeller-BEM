import numpy as np 
import matplotlib.pyplot as plt
from math import *
from Input_Class import Inputs
import toml

class Propeller_Thrust:
    def __init__(self, inp: Inputs):
        self.input = inp 
        self.V = self.input.V #Forward velocity [m/s]
        self.B = self.input.B #Number of blades
        self.R = self.input.R #Radius of the rotor blade [m]
        self.h = self.input.h #Altitude [m]
        self.w = self.input.w #Rotational speed of the rotor [rad/s]
        self.cl = self.input.cl #cl-alpha for alpha = 3deg
        self.cd = self.input.cd #cl-cd for alpha = 3deg
        self.c = self.input.c
        self.a = self.input.a #Induction factor guess
        self.ap = self.input.ap #Axial induction factor guess
        self.solidity = self.input.solidity #Solidity of the blades
        self.n = self.input.n #Number of propellers
        self.graph = self.input.graph

        self.r = None #Radius of blade elements
        self.rho = None #Density at altitude
        self.vi = None #Induced velocity
        self.T_rotor = None
        self.N_rotor = None
        self.D_rotor = None
    
    ############################    ADJUSTING FOR THE ALTITUDE    ############################
    def altitude(self):
        #ISA VALUES
        rho0 = 1.225
        T0 = 273.15 + 15
        a = -0.0065 #K/m
        R = 287 #J/kg*K
        g = 9.80665 #m/s^2

        T1 = T0 + a*self.h
        rho = rho0*(T1/T0)**(-(g/(a*R))-1)
        self.rho = rho
        return rho

    def iteration(self):
        diff_a = 0.1
        diff_ap = 0.1

        while diff_a > 0.05 and diff_ap > 0.05:
            #INFLOW ANGLE
            inflow = np.arctan((self.V*(1+self.a))/(self.w*self.r*(1-self.ap)))


            ftip = (self.B/2)*(self.R - self.r)/(self.r*np.sin(inflow))
            F = (2/np.pi)*np.arccos(np.exp(-ftip)) #tip correction

            a_new = ((4*F*np.sin(inflow)**2)/(self.solidity*(self.cl*np.cos(inflow) - self.cd*np.sin(inflow))) -1)**(-1)
            ap_new = ((4*F*np.sin(inflow)*np.cos(inflow))/(self.solidity*(self.cl*np.sin(inflow) - self.cd*np.cos(inflow)))+1)**(-1)

            diff_a = (a_new - self.a)/self.a
            diff_ap = (ap_new - self.ap)/self.ap
            self.a = a_new
            self.ap = ap_new
        return

    ############################    INDUCTION FACTORS ITERATION   ############################
    def induction_factors(self):
        radius = []
        a = [] #List of optimal induction factors
        ap = [] #List of optimal axial induction factors

        for r in np.arange(1, 6, 0.05):
            self.r = r 
            self.iteration()
            radius.append(r)
            a.append(self.a)
            ap.append(self.ap)

        if self.graph == True:
            fig, ax = plt.subplots(nrows=1, ncols=2)
            ax[0].plot(radius, a, label = 'Axial induction factor')
            ax[0].grid(True)
            ax[0].set_xlabel('Radius of the blade [m]')
            ax[0].set_ylabel('Axial induction factor [m]')

            ax[1].plot(radius, ap, label = 'Tangential induction factor')
            ax[1].grid(True)
            ax[1].set_xlabel('Radius of the blade [m]')
            ax[1].set_ylabel('Tangential induction factor')

            
            ax[0].set_title('a vs r')
            ax[1].set_title(r'a prime vs r')
            
            fig.suptitle(f'Axial & Tangential Induction factors for a Speed of {np.round(self.V*3.6)} km/hr')
            plt.show()

        self.r = radius
        self.a = a
        self.ap = ap
        return

    ############################       THRUST OF THE ROTOR       #############################
    def thrust(self):
        T = []
        Q = []
        D = []
        V_ind = []
        F_cor = []
        twist_opt = []
        for i in np.arange(0, len(self.a)):
            dr = 0.05
            #INFLOW ANGLE
            inflow = np.arctan((self.V*(1+self.a[i]))/(self.w*self.r[i]*(1-self.ap[i])))

            ftip = (self.B/2)*(self.R - self.r[i])/(self.r[i]*np.sin(inflow))
            F = (2/np.pi)*np.arccos(np.exp(-ftip)) #tip correction

            dT = 4*np.pi*(self.r[i])*self.rho*((self.V)**2)*(1+self.a[i])*self.a[i]*F*dr

            dQ = self.rho*4*np.pi*((self.r[i])**3)*(self.V)*(1+self.a[i])*self.ap[i]*self.w*F*dr

            V0 = self.V*(1+self.a[i])
            V2 = self.w*self.r[i]*(1-self.ap[i])

            V1 = np.sqrt(V0**2 + V2**2)
            dD = self.B*self.rho*0.5*(V1**2)*self.cd*0.5
            Vind = self.V*(1+self.a[i]*F)
            
            alpha = 3 #degrees
            inflow = inflow*(180/np.pi)
            twist_opt.append(alpha+inflow)
            F_cor.append(F)
            T.append(dT)
            Q.append(dQ)
            D.append(dD)

        Pr = (self.n*sum(Q)*self.w + self.n*sum(D))/1000 
        Q = np.round(sum(Q)*self.n/1000, 2)
        T = np.round(sum(T)*self.n/1000, 2)

        if self.graph == True: 
            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.plot(self.r, F_cor)
            ax1.set_xlabel('Radius [m]')
            ax1.set_ylabel('Correction factor [-]')
            ax1.set_title('Correction factor vs Radius')
            ax2.plot(self.r, V_ind)
            ax2.set_title('Axial Velocity vs Radius')
            ax2.set_xlabel('Radius [m]')
            ax2.set_ylabel('Axial induced velocity [m/s]')
            fig.suptitle('Propeller Effects')
            ax1.grid(True)
            ax2.grid(True)
            plt.show()

        return Pr, T, Q, twist_opt
    
    def callallnames(self):
        self.altitude()
        self.induction_factors()
        Pr, T, Q, twist_opt = self.thrust()

        return Pr, T, Q, twist_opt
    
if __name__ == '__main__':
   ############################       REPLACE INPUTS HERE       ############################# 
    file = 'Propeller-BEM\\Input_information.toml'
    with open(file, 'r') as f:
        data = toml.load(f)

    p = Propeller_Thrust(Inputs(data))
    Pr, T, Q, twist_opt = p.callallnames()
    
    print(f'A blade twist of {np.round(max(twist_opt)-min(twist_opt))} degrees is optimal.')
    print(f'The total amount of power required is {np.round(Pr, 2)} kW')
    print(f'The total amount of thrust produced is {np.round(T, 2)} kN')

