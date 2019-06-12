# -*- coding: utf-8 -*-

import re
from sys import argv,path
path.insert(0,r'../../')
from JorG.loadsave import POSCARloader
import numpy as np
from itertools import product

class error:
    unexcepted = 12

BRUTAL='Fe'
CUTOFF=2.46
NUMBER=1

class MAGMOMloader:
    rawTxt = []
    data   = []
    def __init__(self,*args,**kwargs):
        for inputName in args:
            try:
                with open(inputName,"r+") as inFile:
                    self.rawTxt.append(inFile.readlines())
            except FileNotFoundError:
                print("File \"%s\" not found!"%inputName)
            except Exception:
                print("Unexcepted error!")
                exit(error.unexcepted)

    def __len__(self):
        return len(self.data)

    @staticmethod
    def should_skip(line):
        for regex in ['^\s+$','^\s*#','^-+$',
                      '^\s*tot','magnetization \(x\)']:
            if re.search(regex,line):
                return True
        return False

    @staticmethod
    def should_stop__reading(line,ISTOBEREAD):
        if re.search('^\s*tot',line):
            return False
        return ISTOBEREAD

    @staticmethod
    def should_start_reading(line,ISTOBEREAD):
        if re.search(' \(x\)',line):
            return True
        return ISTOBEREAD

    @staticmethod
    def read_energy(line):
        if "energy  without entropy" in line:
            return float(line.split()[-1])

    @staticmethod
    def read_moments(text):
        moments={}
        energy = None
        ISTOBEREAD=False
        for line in text:
            ISTOBEREAD = MAGMOMloader.should_start_reading(line,ISTOBEREAD)
            ISTOBEREAD = MAGMOMloader.should_stop__reading(line,ISTOBEREAD)
            if MAGMOMloader.should_skip(line):
                continue
            elif ISTOBEREAD:
                moments[int(line.split()[0])] = float(line.split()[4])
            elif energy is None:
                energy = MAGMOMloader.read_energy(line)
        return {'moments': moments, 'energy': energy}

    def parse_text(self,text):
        self.data.append(MAGMOMloader.read_moments(text))

    def parse(self):
        self.data   = []
        for text in self.rawTxt:
            self.parse_text(text)

    def get_energy(self,idx=0):
        return self.data[idx]['energy']

    def get_moments(self,idx=0):
        return self.data[idx]['moments']

    def __call__(self,idx=0):
        return self.data[idx]

#
#░█▀▀░█▄█░█▀█░█▀▄░▀█▀░░░█▀█░▀█▀░█▀▀░█░█░░░░░█░█░█▀█
#░▀▀█░█░█░█▀█░█▀▄░░█░░░░█▀▀░░█░░█░░░█▀▄░▄▄▄░█░█░█▀▀
#░▀▀▀░▀░▀░▀░▀░▀░▀░░▀░░░░▀░░░▀▀▀░▀▀▀░▀░▀░░░░░▀▀▀░▀░░
#
class SmartPickUp:
    POSCARs           = None
    MAGMOMs           = None
    systemOfEquations = []
    def __init__(self,   numberOfNeighbors,
                 namesOfInteractingAtoms):
        self.numberOfNeighbors       = numberOfNeighbors
        self.namesOfInteractingAtoms = namesOfInteractingAtoms
        self.solution = None

    def read_POSCARs(self,*args):
        self.solution = None
        POSCARs = [ "%s/POSCAR"%arg for arg in args ]
        self.POSCARs = POSCARloader(*POSCARs,spam=False)
        self.POSCARs.parse()

    def read_MAGMOMs(self,*args):
        self.solution = None
        OUTCARs = [ "%s/OUTCAR"%arg for arg in args ]
        self.MAGMOMs = MAGMOMloader(*OUTCARs,spam=False)
        self.MAGMOMs.parse()

    def read(self,*args):
        self.solution = None
        self.read_POSCARs(*args)
        self.read_MAGMOMs(*args)

    def make_crystal(self,idx=0):
        self.solution = None
        self.crystal = []
        for mul,(i,atom) in product(product(range(-1,1),repeat=3),
                                    enumerate(self.POSCARs(idx)['cell'])):
            boost = np.dot(mul,self.POSCARs(idx)['directions'])
            self.crystal.append([i,atom[0],atom[1]+boost])

    def map_distances(self,idx=0):
        self.solution = None
        self.distances = set([])
        self.make_crystal(idx)
        for atom in self.POSCARs(idx)['cell']:
            d = np.around(np.linalg.norm(atom[1]),decimals=2)
            if atom[0] in self.namesOfInteractingAtoms:
                self.distances.add(d)

        self.distances = np.array([d for d in self.distances])
        self.distances = np.sort(self.distances)
        self.distances = self.distances[1:1+self.numberOfNeighbors]

    # for sorted!
    def get_system_of_equations(self):
        self.solution = None
        self.map_distances()

        self.systemOfEquations = []
        E0 = self.MAGMOMs(0)['energy']
        self.dE = []
        for i in range(1,len(self.MAGMOMs)):
            try:
                self.dE.append(self.MAGMOMs(i)['energy']-E0)
            except TypeError:
                print("VASP hasn't finished this run (%d/%d)"%(i,len(self.MAGMOMs)-1))
                continue
            self.systemOfEquations.append(np.zeros(self.numberOfNeighbors))
            flipped = []
            for idx,atom in enumerate(self.POSCARs(0)['cell']):
                if atom[0] not in self.namesOfInteractingAtoms:
                    continue

                momentA = self.MAGMOMs(0)['moments'][idx+1]
                momentB = self.MAGMOMs(i)['moments'][idx+1]
                scalar = momentA * momentB
                if (abs(scalar) > 1e-5 and scalar < 0.0):
                    flipped.append(idx)
        
            for f,(j,atom) in product(flipped,enumerate(self.crystal)):
                if atom[1] not in self.namesOfInteractingAtoms:
                    continue
                elif atom[0] in flipped:
                    continue
                d = np.around(np.linalg.norm(atom[2]-self.POSCARs(0)['cell'][f][1]),decimals=2)
                if d in self.distances:
                    self.systemOfEquations[i-1][np.argwhere(self.distances==d)] -= 2*self.MAGMOMs(0)['moments'][atom[0]+1]*self.MAGMOMs(i)['moments'][f+1]

    energyRatios = {'eV' :    1.0,
                    'meV':    1000.0,
                    'Ry' :    np.reciprocal(13.6056980659),
                    'mRy':    1000.0*np.reciprocal(13.6056980659),
                    'He' :    0.5*np.reciprocal(13.6056980659),
                    'mHe':    500.0*np.reciprocal(13.6056980659),
                    'K'  :    11604.51812}

    @staticmethod
    def convert(*args,**kwargs):
        try:
            return [SmartPickUp.energyRatios[kwargs['units']]*arg for arg in args]
        except KeyError:
            print("No unit defined! Values will be in eV")
            return args

    def solve_system_of_equations(self, ISREDUNDANT=False, **kwargs):
        self.systemOfEquations = np.array(self.systemOfEquations)
        self.dE                = np.array(self.dE)

        if len(self.systemOfEquations) == 1:
            self.solution = self.dE[0]/self.systemOfEquations[0]
        elif self.systemOfEquations.shape[0] == self.systemOfEquations.shape[1]:
            self.solution = np.linalg.solve(self.systemOfEquations,self.dE)
        else:
            self.solution = np.linalg.lstsq(self.systemOfEquations,self.dE,rcond=None)[0]

    def solve(self,**kwargs):
        if self.solution is not None:
            return SmartPickUp.convert(*self.solution,**kwargs)
        self.get_system_of_equations()
        self.solve_system_of_equations(**kwargs)
        return SmartPickUp.convert(*self.solution,**kwargs)
