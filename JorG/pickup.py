# -*- coding: utf-8 -*-

import re
from sys import argv,path
path.insert(0,r'../../')
from JorG.loadsave import POSCARloader
from JorG.generator import apply_mirrorsXYZ
import numpy as np
from itertools import product
from copy import copy

class error:
    unexcepted = 12

class EnergyConverter:
    energyRatios = {'eV' :    1.0,
                    'meV':    1000.0,
                    'Ry' :    np.reciprocal(13.6056980659),
                    'mRy':    1000.0*np.reciprocal(13.6056980659),
                    'He' :    0.5*np.reciprocal(13.6056980659),
                    'mHe':    500.0*np.reciprocal(13.6056980659),
                    'K'  :    11604.51812}
    default = { 'groundMoment'  : 1.0,
                'excitedMoment' : 1.0}
    types = [ "       no moment",
              "  average moment",
              "geometric moment",
              " original moment",
              " neoteric moment" ]   
    @staticmethod
    def convert(*args,**kwargs):
        # Returns array of J values with different conventions (see types)
        settings = EnergyConverter.default
        settings.update(kwargs)
        try:
            notMomentSq = 1.0
            avgMomentSq = 0.25*(settings['groundMoment']+settings['excitedMoment'])**2
            geoMomentSq = settings['groundMoment']*settings['excitedMoment']
            orgMomentSq = settings['groundMoment']*settings['groundMoment']
            newMomentSq = settings['excitedMoment']*settings['excitedMoment'] 
            return [ [notMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args],
                     [avgMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args],
                     [geoMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args],
                     [orgMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args],
                     [newMomentSq * EnergyConverter.energyRatios[settings['units']]*arg for arg in args]]

        except KeyError:
            print("No unit defined! Values will be in eV")
            return args

class ReadMoments:
    @staticmethod
    def should_skip(line):
        for regex in ['^\s+$','^\s*#','^-+$',
                      '^\s*tot','magnetization \(x\)']:
            if re.search(regex,line):
                return True
        return False

    def should_stop__reading(self,line):
        if re.search('^\s*tot',line):
            self.ISTOBEREAD = False

    def should_start_reading(self,line):
        if re.search(' \(x\)',line):
            self.ISTOBEREAD = True

    @staticmethod
    def read_energy(line):
        if "energy  without entropy" in line:
            return float(line.split()[-1])

    def __init__(self):
        self.moments    = {}
        self.energy     = None
        self.ISTOBEREAD = False

    def __call__(self,text):
        for line in text:
            self.read_line(line)
        return {'moments': self.moments, 'energy': self.energy}

    def read_line(self,line):
        self.should_start_reading(line)
        self.should_stop__reading(line)
        if self.should_skip(line):
            return
        if self.ISTOBEREAD:
            self.moments[int(line.split()[0])] = float(line.split()[4])
        elif self.energy is None:
            self.energy = self.read_energy(line)

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

    def parse_text(self,text):
        read_moments = ReadMoments()
        self.data.append(read_moments(text))

    def parse(self):
        self.data   = []
        for text in self.rawTxt:
            self.parse_text(text)

    def get_energy(self,idx=0):
        return self.data[idx]['energy']

    def get_moments(self,idx=0):
        return self.data[idx]['moments']

    def get_average_magnitude(self,idx=0,indices=None):
        if indices is not None:
            return np.average(np.abs(np.take(list(self.data[idx]['moments'].values()),indices)))
        return np.average(np.abs(np.array(list(self.data[idx]['moments'].values()))))

    def __call__(self,idx=0):
        return self.data[idx]

#
#░█▀▀░█▄█░█▀█░█▀▄░▀█▀░░░█▀█░▀█▀░█▀▀░█░█░░░░░█░█░█▀█
#░▀▀█░█░█░█▀█░█▀▄░░█░░░░█▀▀░░█░░█░░░█▀▄░▄▄▄░█░█░█▀▀
#░▀▀▀░▀░▀░▀░▀░▀░▀░░▀░░░░▀░░░▀▀▀░▀▀▀░▀░▀░░░░░▀▀▀░▀░░
#
class SmartPickUp:
    def __init__(self,   numberOfNeighbors,
                 namesOfInteractingAtoms):
        self.numberOfNeighbors       = numberOfNeighbors
        self.namesOfInteractingAtoms = namesOfInteractingAtoms
        self.types                   = EnergyConverter.types

    def read_POSCARs(self,*args):
        lastBackSlashRemoved = [ re.sub('/$','',arg) for arg in args ]
        POSCARs = [ "%s/POSCAR"%arg for arg in lastBackSlashRemoved ]
        self.POSCARs = POSCARloader(*POSCARs,spam=False)
        self.POSCARs.parse()

    def read_MAGMOMs(self,*args):
        OUTCARs = [ "%s/OUTCAR"%arg for arg in args ]
        self.MAGMOMs = MAGMOMloader(*OUTCARs,spam=False)
        self.MAGMOMs.parse()

    def read(self,*args,**kwargs):
        self.read_MAGMOMs(*args)
        self.read_POSCARs(*args)
        if 'reference' in kwargs:
            self.reference = self.POSCARs(0)['cell'][kwargs['reference']][1]
            self.ref       = kwargs['reference']
        else:
            print("Warning: reference @ 0. Is that ok?")
            self.ref       = 0
            self.reference = self.POSCARs(0)['cell'][0][1]

    def make_crystal(self,idx=0):
        self.crystal  = self.POSCARs(idx)['cell']
        self.crystal  = [ [atom[0],atom[1],self.MAGMOMs.get_moments()[i+1]] for i,atom in enumerate(self.crystal) ]
        self.crystal8 = apply_mirrorsXYZ(self.POSCARs(0)['directions'],self.crystal,
                                         reference=self.ref)

    def map_distances(self,idx=0):
        self.distances = set([])
        self.make_crystal(idx)
        for atom in self.POSCARs(idx)['cell']:
            d = np.around(np.linalg.norm(atom[1]-self.reference),decimals=2)
            if atom[0] in self.namesOfInteractingAtoms:
                self.distances.add(d)
        self.distances = np.sort(np.array(list(self.distances)))[1:1+self.numberOfNeighbors]

    # for sorted!
    def get_system_of_equations(self):
        self.map_distances()
        self.systemOfEquations      = []
        self.dE                     = []
        self.flippingConfigurations = []
        for i in range(1,len(self.MAGMOMs)):
            try:
                self.dE.append(self.MAGMOMs(i)['energy']-self.MAGMOMs(0)['energy'])
            except TypeError:
                print("VASP hasn't finished this run (%d/%d)"%(i,len(self.MAGMOMs)-1))
                continue
            self.set_flipps(i)
        self.model             = NaiveHeisenberg(self.flippingConfigurations,self.crystal,self.crystal8)
        self.model.MAGMOMs     = self.MAGMOMs
        self.flipped           = np.unique(np.where(self.flippingConfigurations)[1])
        self.systemOfEquations = self.model.generate(self.namesOfInteractingAtoms,self.distances)
        self.solver            = EquationSolver(self.systemOfEquations,self.dE)
        self.solver.remove_tautologies()

    def set_flipps(self,i):
        self.flippingConfigurations.append([])
        for idx,atom in enumerate(self.POSCARs(0)['cell']):
            self.get_flip(i,idx,atom)

    def get_flip(self,i,idx,atom):
        if atom[0] not in self.namesOfInteractingAtoms:
           self.flippingConfigurations[-1].append(False)
           return
        momentA = self.MAGMOMs(0)['moments'][idx+1]
        momentB = self.MAGMOMs(i)['moments'][idx+1]
        scalar = momentA * momentB
        if (abs(scalar) > 1e-5 and scalar < 0.0):
            self.flippingConfigurations[-1].append(True)
            return
        self.flippingConfigurations[-1].append(False)

    def solve(self,**kwargs):
        try:
            self.solver
        except AttributeError:
            self.get_system_of_equations()
        averageMground = self.MAGMOMs.get_average_magnitude(0,self.flipped)
        averageMexcite = np.average([self.MAGMOMs.get_average_magnitude(i,self.flipped) for i in range(1,len(self.POSCARs))])
        return EnergyConverter.convert(*(self.solver.solve()),groundMoment=averageMground,excitedMoment=averageMexcite,**kwargs)

class EquationSolver:
    def __init__(self,equations,vector):
        self.equations = np.array(equations)
        self.vector    = np.array(vector)
        self.solution  = None

    def solve(self,**kwargs):
        if self.solution is not None:
            return self.solution
        self.solve_system_of_equations(**kwargs)
        return self.solution

    def solve_system_of_equations(self, **kwargs):
        if len(self.equations) == 1:
            self.solution = self.vector[0]/self.equations[0]
        elif self.equations.shape[0] == self.equations.shape[1]:
            self.solution = np.linalg.solve(self.equations,self.vector,**kwargs)
        else:
            self.solution = np.linalg.lstsq(self.equations,self.vector,rcond=None,**kwargs)[0]

    def remove_tautologies(self, **kwargs):
        # removing 0 = 0 equations !
        scale       = np.sum(
                        np.abs(self.equations))
        tautologies = np.argwhere(
                        np.apply_along_axis(
                          np.linalg.norm,1,self.equations)/scale<1e-3)[:,0]
        self.equations = np.delete(self.equations,tuple(tautologies),axis=0)
        return self.equations

    @staticmethod
    def is_to_be_removed(i,j,equations):
        if i == j:
            return False
        scale = np.sum(np.abs(equations))
        inner_product = np.inner(equations[i],equations[j])
        norm_i = np.linalg.norm(equations[i])
        norm_j = np.linalg.norm(equations[j])
        if np.abs(inner_product - norm_j * norm_i)/scale < 1E-8:   
            return True
        return False

    def remove_linears(self, **kwargs):
        # Based on https://stackoverflow.com/questions/
        #                  28816627/how-to-find-linearly-independent-rows-from-a-matrix
        # We remove lineary dependent rows
        remover=set([])
        for i,j in product(range(self.equations.shape[0]),repeat=2):
            if i in remover:
                continue
            if self.is_to_be_removed(i,j,self.equations):
                remover.add(j)
        if remover:
            self.equations = np.delete(self.equations,tuple(remover),axis=0)
            return remover

    def remove_linear_combinations(self, secondSystem):
        scale            = np.sum(np.abs(self.equations))
        resultantSystem  = np.array([self.equations[0]])
        gram_schmidt     = np.array([self.equations[0]])
        self.equations   = np.delete(self.equations, (0), axis=0)
        resultantSystem2 = np.array([secondSystem[0]])
        secondSystem     = np.delete(secondSystem, (0), axis=0)
        while len(resultantSystem) < self.equations.shape[1]:
            if self.equations.size == 0:
                print("ERROR! Not enough equations! Please rerun.")
                exit(-3)
            tmpVector = np.copy(self.equations[0])
            for vector in gram_schmidt:
                tmpVector -= np.inner(tmpVector,vector)*vector/np.inner(vector,vector)
            if np.linalg.norm(tmpVector)/scale > 1e-4:
                gram_schmidt     = np.vstack((gram_schmidt,tmpVector))
                resultantSystem  = np.vstack((resultantSystem,self.equations[0]))
                self.equations   = np.delete(self.equations, (0), axis=0)
                resultantSystem2 = np.vstack((resultantSystem2,secondSystem[0]))
                secondSystem     = np.delete(secondSystem, (0), axis=0)
        self.equations = resultantSystem
        return resultantSystem,resultantSystem2

class DefaultMoments:
    def __init__(self,crystal):
        self.moments = [ atom[2] for atom in crystal ]
        self.moments.insert(0,0.0)
    def __call__(self,idx=0):
        return {'energy': 0.0, 'moments': self.moments}

class NaiveHeisenberg:
    def __init__(self,flippings,
                 crystal,crystal8):
        self.MAGMOMs   = DefaultMoments(crystal)
        self.flippings = flippings
        self.crystal   = crystal
        self.crystal8  = crystal8

    def generate(self,mask,flipper):
        self.flipper = np.array(flipper)
        self.mask    = mask
        self.systemOfEquations = np.zeros((len(self.flippings),len(flipper)))
        for (i,config),(I,atomI),atomJ in product(enumerate(self.flippings),
                                                  enumerate(self.crystal),
                                                            self.crystal8):
            if config[I] == config[atomJ[3]]:
                continue
            distance = self.check_if_contributes(config,atomI,atomJ)
            if not distance: 
                continue
            j = np.argwhere(np.abs(self.flipper - distance)<1e-2)
            if j.size:
                self.systemOfEquations[i][j[0][0]] += 2.0*atomI[2]*self.MAGMOMs(i)['moments'][atomJ[3]+1]
        return self.systemOfEquations

    def check_if_contributes(self,config,atomI,atomJ):
        if (atomI[2] != 0.0 and atomI[0] in self.mask
        and atomJ[2] != 0.0 and atomJ[0] in self.mask):
            distance = np.linalg.norm(atomI[1]-atomJ[1])
            if np.abs(distance) > 1e-2:
                return distance
        return False

class Reference:
    def __init__(self,POSCAR):
        loader = POSCARloader(POSCAR)
        loader.parse()
        firstLine = loader()['comment']
        try:
            self.reference = int(re.search('NewRef:\s*([0-9]+),',firstLine).group(1))
        except AttributeError as err:
            print("Reference not found in POSCAR comment! Taking 0!")
            print(err)
            self.reference = 0
        except ValueError as err:
            print("Reference cannot be converted to int! Taking 0!")
            print(err)
            self.reference = 0
        if self.reference < 0:
            print("Wrong reference (%d < 0)! Taking 0!"%self.reference)
            self.reference = 0
    def __str__(self):
        return str(self.reference)
    def __call__(self):
        return self.reference
