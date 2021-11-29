# -*- coding: utf-8 -*-

import numpy as np
import errno
from os.path import exists
from datetime import datetime
from os import makedirs

class Errors:
    failed_to_generate = 201
    no_reference       = 22
    systemerror        = 255
    erroneous_input    = 44

import JorGpi.generate.loadsave as loadsave
from JorGpi.POSCARloader import POSCARloader

class StreamHandler:
    def __init__(self,*args):
        self.streams = []
        for arg in args:
            name = StreamHandler.fix_name(arg)
            StreamHandler.make_directory(name)
            self.streams.append(name)

    @staticmethod
    def fix_name(name):
        if name is None:
            # if output directory is not given:
            name = "output/"+datetime.now().strftime("%Y%m%d%H%M%S")
        else:
            # remove multiple '/' and possible '/' at the end
            name = re.sub('/+','/',name)
            name = re.sub('/$','',name)
        return name

    @staticmethod
    def make_directory(name):
        # creating output path
        try:
            makedirs(name)
        except OSError as err:
            if err.errno == errno.EEXIST:
                print("%s exists: Data will be overwritten!"%name)

    def __call__(self, idx=0):
        return self.streams[idx]

    @staticmethod
    def load_vasp(poscarfile,incarfile):
        if not exists(poscarfile):
            print('POSCAR file %s does not exist!'%poscarfile)
            exit(Errors.erroneous_input)
        if incarfile is not None and not exists(incarfile):
            print('INCAR file %s does not exist!'%incarfile)
            exit(Errors.erroneous_input)

        load_POSCAR          = POSCARloader(poscarfile)
        load_POSCAR.parse()
        readData             = load_POSCAR(0)

        if incarfile is None:
            oldMoments=[ 1.0 ] * len(readData['cell'])
            return readData,oldMoments,''

        load_INCAR           = loadsave.INCARloader(readData['cell'],fileName=incarfile,
                                                    atomNames=readData['atomNames'])
        oldMoments,incarData = load_INCAR()
        return readData,oldMoments,incarData

from copy import copy
class Standard:
    values = {'scale'     : 5.0,
              'vector'    : 0.1,
              'background': [253,246,227],
              'radius'    : 10.0,
              'center'    : [0.0, 0.0, 0.0],
              'color'     : [176,255,176]}
    @staticmethod
    def fix(**kwargs):
        result = copy(Standard.values)
        result.update(kwargs)
        return result

class JmolVisualization:
    scriptText =\
"""vibration ON
vector ON
vector SCALE %.2f
vector %.2f
set bondRadiusMilliAngstroms 0
background [%d,%d,%d]
isoSurface SURF center {%.5f %.5f %.5f} sphere %.5f
color $SURF translucent [%d,%d,%d]"""
    @staticmethod
    def create_script(outDirName,**kwargs):
        kwargs = Standard.fix(**kwargs)
        with open(outDirName+"/script.jmol","w+") as stream:
            stream.write(JmolVisualization.scriptText%(
                kwargs['scale'],
                kwargs['vector'],
                *kwargs['background'],
                *kwargs['center'],
                kwargs['radius'],
                *kwargs['color'],))

from os import remove
class TemporaryFiles:
    names     = [".input",".supercell",".directions"]
    extension = ".dat"
    def __init__(self, suffix=str(np.random.randint(1000,9999)), prefix=''):
        self.suffix = suffix
        self.prefix = prefix

    def write_input(self,allFlippable,crystal):
        name = self.prefix+self.names[0]+self.suffix+self.extension
        with open(name,"w+") as isingFile:
            for i in allFlippable:
                isingFile.write("%d %.8f %.8f %.8f %.2f\n"%(
                i,*crystal[i][1],crystal[i][2]))

    def write_supercell(self,crystal):
        name = self.prefix+self.names[1]+self.suffix+self.extension
        with open(name,"w+") as supercellFile:
            for i,atom in enumerate(crystal):
                supercellFile.write("%d %.8f %.8f %.8f %.2f\n"%(
                i,*atom[1],atom[2]))

    def write_directions(self,directions):
        name = self.prefix+self.names[2]+self.suffix+self.extension
        with open(name,"w+") as dirFile:
            for direction in directions:
                dirFile.write("%.8f %.8f %.8f\n"%(*direction,))

    def __str__(self):
        return "%s %s %s"%(self.get_files())

    def get_files(self):
        iName = self.prefix+self.names[0]+self.suffix+self.extension
        sName = self.prefix+self.names[1]+self.suffix+self.extension
        dName = self.prefix+self.names[2]+self.suffix+self.extension
        return dName,sName,iName

    def __del__(self):
        for name in self.names:
            try:
                remove(self.prefix+name+self.suffix+self.extension)
            except OSError:
                print("No file %s"%(self.prefix+name+self.suffix+self.extension))

from JorGpi.aux.PeriodicTable import periodicTableElement
import re
import spglib
class VariableFixer:
    @staticmethod
    def fix_reference(reference,cell,mask):
        if reference >= 0:
            return cell[reference],reference
        for i,atom in enumerate(cell):
            if "$"+atom[0]+"$" in mask:
                return atom,i
        print("Error:\n%s\nare not in input file!"%re.sub('\$',' ',mask))
        exit(Errors.no_reference)

    @staticmethod
    def from_refined(refinedCell):
        directions = np.array(refinedCell[0])
        cell = []
        for refinedAtom,atomType in zip(refinedCell[1],refinedCell[2]):
            newPosition = np.dot(refinedAtom,refinedCell[0])
            cell.append((periodicTableElement[atomType-1],newPosition))
        return cell,directions

    @staticmethod
    def fix_neighbor(nearestNeighbor,cutOff):
        # by defauld: second NN
        if nearestNeighbor is None:
            if cutOff is None:
                return 2,None
            return -1,cutOff
        return nearestNeighbor,None

    @staticmethod
    def fix_directions(copiesInEachDirection,directions):
        return [(mul+1)*d for mul,d in zip(copiesInEachDirection, directions)]

    @staticmethod
    def add_to_all(arr,add=1):
        return [ a+add for a in arr]

class Symmetry:
    def __init__(self,cellSymmetry):
        self.cellSymmetry = cellSymmetry
        self.symmetry     = self.get_symmetry(self.cellSymmetry)

    @staticmethod
    def get_symmetry(cell):
        return spglib.get_symmetry_dataset(cell)

    def standarize(self):
        return (spglib.standardize_cell(self.cellSymmetry, to_primitive=1,
                                        no_idealize=0, symprec=1e-1))

    def get_standarized(self):
        refinedCell = (spglib.refine_cell(self.cellSymmetry,
                                              symprec=1e-1))
        return spglib.get_symmetry_dataset(self.standarize()),\
               spglib.get_symmetry_dataset(refinedCell)

from JorGpi.aux.format import Color, print_vector, print_label
from JorGpi.aux.format import print_crystal, print_moments
class Msg:
    @staticmethod
    def print_equations(equations,isNotRedundant=False):
        print_label("System of equations:",labelStyle=Color.bold)
        for equation in equations:
            print_vector(equation)
        if not isNotRedundant:
            print_label("Redundant system of equations.",
                        labelStyle=Color.bold)
            print_label("Least square method is to be used to obtain Heisenberg model.",
                        labelStyle=Color.bold)
            print_label("It may be better. But it may also mess everything.",
                        labelStyle=Color.bold)
        else:
            print_label("det SoE = %.1e"%np.linalg.det(equations),
                        labelStyle=Color.bold)

    @staticmethod
    def print_solver_status(configs):
        if configs > 1e6:
            try:
                print_label("Checking total number of configurations: %.2e"%float(configs),
                                labelStyle=Color.bold+Color.darkred)
            except OverflowError:
                print_label("Checking total number of configurations: incomprehensible",
                                labelStyle=Color.bold+Color.darkred)
        else:
            print_label("Checking total number of configurations: %d"%configs,
                            labelStyle=Color.bold+Color.darkred)
        print_label("Preparing and running ASA-solver...",
                      labelStyle=Color.bold+Color.blue)

    @staticmethod
    def print_crystal_info(**kwargs):
        try:
            print_label(kwargs['title'])
        except KeyError:
            print("placeholder",end="")
        try:
            print_label("Size: %dx%dx%d"%(kwargs['copies']),labelStyle=Color.bold)
        except KeyError:
            print_label("Size: 1x1x1",labelStyle=Color.bold)
        try:
            print_crystal(kwargs['crystal'],directions=kwargs['directions'])
        except KeyError:
            print("",end="")
        try:
            print_label("Reference atom in the system is No. %d:"%(kwargs['reference']+1),
                        atoms=[kwargs['crystal'][kwargs['reference']]],
                        vectorStyle=Color.darkcyan,labelStyle=Color.bold)
        except KeyError:
            print("",end="")
        try:
            print_moments(kwargs['moments'],cell=kwargs['crystal'])
        except KeyError:
            print("",end="")

def read_flips(name='best.flips'):
    try:
        read = np.loadtxt(name,bool)
        return read
    except OSError:
        print("No %s here..."%name)
        return
