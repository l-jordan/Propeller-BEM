import toml
import numpy as np

class Inputs():
    def __init__(self, data) -> None:
        super().__init__()

        self.graph = data['Graph']['graph']
    
        self.V = data['Aircraft Environment']['V'] #velocity [m/s]
        self.h = data['Aircraft Environment']['h'] #altitude [m]

        self.B = data['Propeller information']['B']
        self.R = data['Propeller information']['R']
        self.rh = data['Propeller information']['rh']
        self.rt = data['Propeller information']['rt']
        self.RPM = data['Propeller information']['RPM']
        self.w = self.RPM*(2*np.pi/60)
        self.c = data['Propeller information']['chord']
        self.n = data['Propeller information']['Prop number']

        self.cl = data['Airfoil']['cl']
        self.cd = data['Airfoil']['cd']

        self.a = data['Induction factors']['a']
        self.ap = data['Induction factors']['ap']

        self.rm = ((self.rt**2 + self.rh**2)/2)**(0.5)
        self.s = (2*np.pi*self.rm)/self.B
        self.solidity = self.c/self.s


