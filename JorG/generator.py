# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import numpy as np
def get_number_of_pictures(directions,cutOff,referenceAtom=[0,np.zeros(3)]):
    """
        Finding the amount of copies
        of the original SuperCell
                                        """
    multipliers = [] # returned 
    dDirs = np.tile(directions,(2,1)) 
    #                                                                  
    #          .---------------------------.    
    #         / \                           \   
    #        /   \                           \  
    #       /     \                           \ 
    #      /       ^---------------------------.
    #     /       /                           /  
    #    /       /     reference point       /    
    #   /       / C   .                     /     
    #  ^       /                           /       
    # B \     /                           /        
    #    \   /                           /         
    #     \ /                           /          
    #      0--------------------------->           
    #                   A                                     
    # A,B,C -> Vectors
    # normalA = (BxC)/|BxC|
    # normalB = (CxA)/|CxA|
    # normalC = (AxB)/|AxB|

    for i in range(3):
        normal = np.cross(dDirs[i+1],dDirs[i+2])
        normal /= np.linalg.norm(normal)
        height = np.dot(dDirs[i],normal)
        relative = np.dot(dDirs[i]-referenceAtom[1],normal)
        if cutOff > relative:
            multipliers.append(int((cutOff-relative)/height)) # calculating multipliers
        else:
            multipliers.append(0)
    
    return multipliers

#
#
#
#
#

import numpy as np
from aux.periodic import *
from itertools import product
def generate_crystal(multiplyers,cell,directions,atomNames,reference,moments=None):
    """
        Generator of minimal required SuperCell
                                                """
    try:
        if len(moments) != len(cell):
            moments = None
    except:
        pass

    crystal = []
    for name in atomNames:
        for (x,y,z) in product(range(multipliers[0]+1),
                               range(multipliers[1]+1),
                               range(multipliers[2]+1)):
            if moments is None:
                for atom in cell:
                    if atomNames[atom[0]] == name:
                        position = np.copy(atom[1])
                        for a,n in zip([x,y,z],directions):
                            position += a*n
                            
                        flag = "%s"%atomNames[int(atom[0])]
                        crystal.append([flag,position,elementMagneticMoment[periodicTableElement[atomNames[int(atom[0])]]]])    
            else:
                for atom,moment in zip(cell,moments):
                    if atomNames[atom[0]] == name:
                        position = np.copy(atom[1])
                        for a,n in zip([x,y,z],directions):
                            position += a*n
                            
                        flag = "%s"%atomNames[int(atom[0])]
                        crystal.append([flag,position,moment])

    newReference = None
    for i,atom in enumerate(crystal):
      if np.linalg.norm(atom[1] - cell[reference][1]) < 1e-2:
        newReference = i

    return (crystal,newReference)

#
#
#
#
#

import numpy as np
import spglib
from aux.periodic import periodicTableNumber
def wyckoffs_dict(originalCell,      cell,
                  originalDirections,directions,atomNames):
    symmetryCell = None
    try:
      symmetryCell = (directions,
                      [np.dot(row[1],np.linalg.inv(directions))
                          for row in cell],
                      [periodicTableNumber[atomNames[row[0]]]
                          for row in cell])
    except: 
      symmetryCell = (directions,
                      [np.dot(row[1],np.linalg.inv(directions))
                          for row in cell],
                      [periodicTableNumber[row[0]]
                          for row in cell])

    symmetry = spglib.get_symmetry_dataset(symmetryCell)

    originalSymCell = None
    try:
      originalSymCell = (originalDirections,
                         [np.dot(row[1],np.linalg.inv(originalDirections))
                             for row in originalCell],
                         [periodicTableNumber[atomNames[row[0]]]
                             for row in originalCell])
    except:
      originalSymCell = (originalDirections,
                         [np.dot(row[1],np.linalg.inv(originalDirections))
                             for row in originalCell],
                         [periodicTableNumber[row[0]]
                             for row in originalCell])
    originalSym = spglib.get_symmetry_dataset(originalSymCell)

    wyckoffPositionDict = {}
    for i,atom in enumerate(cell):
        for j,originalAtom in enumerate(originalCell):
            if np.linalg.norm(atom[1]-originalAtom[1]) < 1e-3:
                wyckoffPositionDict[symmetry['wyckoffs'][i]] = originalSym['wyckoffs'][j]

    return wyckoffPositionDict,symmetry,originalSym
#
#
#
#
#

import numpy as np
import spglib
from aux.periodic import *
from itertools import product

class generate_from_NN:
    Wyckoffs='abcdefghijklmnopqrstuvwxyz'
    originalSymmetry    = None
    atomTypeMask        = maskFull
    moments             = None 
    extraMultiplier     = np.zeros(3,dtype=int)
    wyckoffPositionDict = None
    distances           = []
    cutOff              = -1
    crystal             = None 
    crystalSymmetry     = None 
    newReference        = None 
    multipliers         = None 
    ISFOUND             = False


    def __init__(self,cell,
                 referenceAtom,
                 directions,
                 nearestNeighbor,
                 atomNames):
        self.cell            = cell
        self.referenceAtom   = referenceAtom
        self.directions      = directions
        self.nearestNeighbor = nearestNeighbor
        self.atomNames       = atomNames

   
        try:
            originalSymmetryCell = (directions,
                                [np.dot(row[1],np.linalg.inv(directions)) for row in cell],
                                [periodicTableNumber[atomNames[row[0]]] for row in cell])
            originalSymmetry = spglib.get_symmetry_dataset(originalSymmetryCell)
        except:
            print("Failed to generate symmetry!")
            pass
   
        diagonal = np.sqrt(np.sum([np.dot(d-referenceAtom[1],d-referenceAtom[1]) for d in directions]))

    def check_in_cell(self,cell,referenceAtom,directions):
        (self.wyckoffPositionDict,
         symmetry,
         self.originalSymmetry) = wyckoffs_dict(self.cell,
                                           cell,
                                           self.directions,
                                           directions,
                                           self.atomNames)
        self.distances = []

    
        for i,(atom,wyck) in enumerate(
                             zip(cell,
                                 symmetry['wyckoffs'])):
            wyckoff = self.wyckoffPositionDict[wyck]
            atomType = self.atomNames[atom[0]]
            distance = np.around(np.linalg.norm(atom[1] - referenceAtom[1]),2)
            if ( distance     not in self.distances    and
                 wyckoff          in self.Wyckoffs     and
                 "$"+atomType+"$" in self.atomTypeMask    ):
                self.distances.append(distance)
    
        self.distances.sort()
        if np.ma.masked_greater(self.distances,self.cutOff).count() <= self.nearestNeighbor:
            return False
        return True

    def __call__(self):
        if self.ISFOUND:
            return (self.cutOff,
                    self.crystal,
                    self.crystalSymmetry,
                    self.newReference,
                    self.multipliers,
                    self.wyckoffPositionDict)
        self.ISFOUND= self.check_in_cell(self.cell, self.referenceAtom, self.directions)

    
        self.crystal = []
        if len(self.distances) > 1:
            minDirection = self.distances[-1]
        else:
            minDirection = 0.99*np.min([np.linalg.norm(d) for d in self.directions])
        for cutOff in minDirection*np.sqrt(np.arange(1.0,np.power(self.nearestNeighbor+1,3),1.0)):
            self.cutOff = cutOff
            self.multipliers = get_number_of_pictures(self.directions,
                                                 cutOff,
                                                 self.referenceAtom)
            self.multipliers += self.extraMultiplier
            extraDirections = [(mul+1)*d 
                               for mul,d in
                               zip(self.multipliers,
                                   self.directions)]
            self.crystal = []                   
            for name in self.atomNames:
                for (x,y,z) in product(range(self.multipliers[0]+1),
                                       range(self.multipliers[1]+1),
                                       range(self.multipliers[2]+1)):
                               
                    if self.moments is None:
                         for atom in cell:
                             if atomNames[atom[0]] == name:
                                 position = np.copy(atom[1])
                                 for a,n in zip([x,y,z],self.directions):
                                     position += a*n
                                 if isinstance(atom[0],int):
                                     self.crystal.append([atom[0],position,elementMagneticMoment[periodicTableElement[self.atomNames[atom[0]]]]])    
                                 else:
                                     self.crystal.append([atom[0],position,elementMagneticMoment[periodicTableElement[atom[0]]]])    
                    else:
                        for atom,moment in zip(self.cell,self.moments):
                            if self.atomNames[atom[0]] == name:
                                position = np.copy(atom[1])
                                for a,n in zip([x,y,z],self.directions):
                                    position += a*n
                                self.crystal.append([atom[0],position,moment])    
    
            self.newReference = None
            self.newReferenceAtom = None
            for i,atom in enumerate(self.crystal):
                if np.linalg.norm(atom[1] - self.referenceAtom[1]) < 1e-3:
                    self.newReference = i
                    self.newReferenceAtom = atom
    
            self.ISFOUND = self.check_in_cell(
                           self.crystal,
                           self.newReferenceAtom,
                           extraDirections)
            if(self.ISFOUND):
                symmetryCell = (extraDirections,
                                [np.dot(row[1],np.linalg.inv(extraDirections))
                                    for row in self.crystal],
                                [periodicTableNumber[self.atomNames[row[0]]]
                                    for row in self.crystal])
                self.crystalSymmetry = spglib.get_symmetry_dataset(symmetryCell)
                self.cutOff = self.distances[self.nearestNeighbor]
                for atom in self.crystal:
                    atom[0] = self.atomNames[atom[0]]
                return (self.cutOff,
                        self.crystal,
                        self.crystalSymmetry,
                        self.newReference,
                        self.multipliers,
                        self.wyckoffPositionDict)

        return None

#
#
#
#
#

from itertools import product
def apply_mirrorsXYZ(dimensions,cell, cutOff=-1.0, reference=0):
    outputCell = []
    for p in product([0,-1],repeat=3):
        projection = np.array([p])
        translation = np.dot(projection,dimensions)[0]
        for i,atom in enumerate(cell): 
#            if cutOff > 0:
#                if np.linalg.norm(atom[1]+translation - cell[reference][1]) > cutOff:
#                    continue
            outputCell.append([atom[0],atom[1]+translation,atom[2],i])
    return outputCell
