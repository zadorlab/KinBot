import numpy as np
from ase.calculators.calculator import Calculator, all_changes

class CostFunction(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, coords=None, restart=None, label='costfunction', atoms=None, 
                 tnsr=False, **kwargs):
        Calculator.__init__(self, restart=restart, label=label, atoms=atoms, 
                            tnsr=tnsr, **kwargs)
        self.coords = coords

    def calculate(self, atoms=None, properties=['energy', 'forces'], 
                  system_changes=all_changes, args=None):
        
        x = np.reshape(atoms.positions, 3 * len(atoms.positions))
        self.results['energy'] = self.eval(x)
        self.results['forces'] = np.reshape(self.gradient(x) * -1, (len(atoms.positions), 3 ))
        
    def eval(self, x):
        """
        x is a vector of length 3N with N the number of atoms
        containing the cartesian coordinates [x1, y1, z1, x2, ..., xN, yN, zN]
        """
        e = 0
        for coord in self.coords:
            i = coord[0]
            j = coord[1]
            d = coord[2]
            weight = coord[3]
            dc = (x[3 * i] - x[3 * j]) ** 2 + (x[3 * i + 1] - x[3 * j + 1]) ** 2 + (x[3 * i + 2] - x[3 * j + 2]) ** 2
            if len(coord) == 5:
                if dc < d:
                    # add the one-sided potential
                    e += ((dc - d) * weight) ** 2
            else:
                e += ((dc - d) * weight) ** 2
        return e

    def gradient(self, x):
        grad = np.zeros(len(x))
        for coord in self.coords:
            i = coord[0] #atom1
            j = coord[1] #atom2
            d = coord[2] #value^2
            weight = coord[3]
            dc = (x[3 * i] - x[3 * j]) ** 2 + (x[3 * i + 1] - x[3 * j + 1]) ** 2 + (x[3 * i + 2] - x[3 * j + 2]) ** 2 #current dist between i and j
            if len(coord) == 5:
                if dc < d:
                    grad[3 * i] += 2 * ((dc - d) * weight) * 2 * (x[3 * i] - x[3 * j])
                    grad[3 * i + 1] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 1] - x[3 * j + 1])
                    grad[3 * i + 2] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 2] - x[3 * j + 2])

                    grad[3 * j] += 2 * ((dc - d) * weight) * 2 * (x[3 * i] - x[3 * j]) * -1
                    grad[3 * j + 1] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 1] - x[3 * j + 1]) * -1
                    grad[3 * j + 2] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 2] - x[3 * j + 2]) * -1
            else:
                grad[3 * i] += 2 * ((dc - d) * weight) * 2 * (x[3 * i] - x[3 * j])
                grad[3 * i + 1] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 1] - x[3 * j + 1])
                grad[3 * i + 2] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 2] - x[3 * j + 2])

                grad[3 * j] += 2 * ((dc - d) * weight) * 2 * (x[3 * i] - x[3 * j]) * -1
                grad[3 * j + 1] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 1] - x[3 * j + 1]) * -1
                grad[3 * j + 2] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 2] - x[3 * j + 2]) * -1

        return grad
    
    def hessian(self, x):
        
        hess = np.zeros((len(x),len(x)))
        grad = self.gradient(x)
        for x1 in range(len(x)):
            for x2 in range(len(x)):
                hess[x1][x2] = grad[x1] * grad[x2]

        return hess
