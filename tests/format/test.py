#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
import time
import numpy as np

from aux.format import standard,color
from aux.format import print_vector
from aux.format import print_atom
from aux.format import print_case
from aux.format import print_moments

def main(**args):
    pass

def line():
    print(92*'+')

if __name__ == '__main__':
    tracker  = -(time.time())

    print("Test of kwargs fixing:")
    testkwargs = {'linewidth' : 99, 'nonexistent' : 'yetIexist'}
    print(testkwargs)
    testkwargs = standard.fix(**testkwargs)
    print(testkwargs)
    line()

    print("Test of print_vector(",end="")
    vector = [0.1, 0.11, 0.111]
    print(vector,end=")\n")
    print("\tplain:")
    print_vector(vector)
    print("\twith vectorStyle=color.BF:")
    print_vector(vector,vectorStyle=color.BF)
    print("\twith linewidth=23:")
    print_vector(vector,linewidth=23)
    print("\twith end=\'\\nSTOP\\n\':")
    print_vector(vector,end='\nSTOP\n')
    line()

    print("Test of print_atom(",end="")
    atom = ["Cu",np.zeros(3)]
    print(atom,end=")\n")
    print("\tplain:")
    print_atom(atom)
    print("\twith vectorStyle=color.BF:")
    print_atom(atom,vectorStyle=color.BF)
    print("\twith elementStyle=color.DARKRED")
    print_atom(atom,elementStyle=color.DARKRED)
    print("\twith end=\'\\nSTOP\\n\':")
    print_atom(atom,end='\nSTOP\n')
    line()

    print("Test of print_case(",end="")
    case = ["Cu",np.zeros(3)]
    settings = {'atomID': 12, 'caseID': 1}
    print(case,', **',settings,end=")\n",sep="")
    print("\tplain:")
    print_case(case,**settings)
    print("\twith caseStyle=color.BF:")
    print_case(case,**settings,caseStyle=color.BF)
    print("\twith numberStyle=color.BF:")
    print_case(case,**settings,numberStyle=color.BF)
    print("\twith distanceStyle=color.BF:")
    print_case(case,**settings,distanceStyle=color.BF)
    print("\twith wyckoffPosition=\'a\':")
    print_case(case,**settings,wyckoffPosition='a')
    print("\twith distance=1.23")
    print_case(case,**settings,distance=1.23)
    print("\twith end=\'\\nSTOP\\n\':")
    print_case(case,**settings,end='\nSTOP\n')
    line()

    print("\tTest print_moments(",end="")
    moments = [0.011,0.880,0.011,0.011,1.900,
               4.800,0.011,0.011,0.011,0.011,
               0.011,0.011,0.011]
    cell    = [('C',  np.array([0.0, 0.0, 0.0])),
               ('Ce', np.array([1.0, 0.0, 0.0])),
               ('Cs', np.array([0.0, 1.0, 0.0])),
               ('In', np.array([0.0, 0.0, 1.0])),
               ('Co', np.array([0.0, 1.0, 1.0])),
               ('Cr', np.array([1.0, 0.0, 1.0])),
               ('W',  np.array([1.0, 1.0, 0.0])),
               ('Re', np.array([1.0, 1.0, 1.0])),
               ('K',  np.array([0.5, 0.0, 0.0])),
               ('V',  np.array([0.0, 0.5, 0.0])),
               ('B',  np.array([0.0, 0.0, 0.5])),
               ('Te', np.array([0.5, 0.5, 0.5])),
               ('Cu', np.array([0.7, 0.7, 0.7]))]
    print(moments,", cell =",cell,")",sep="")
    print("\twit     cell:")
    print_moments(moments,cell=cell)
    print("\twithout cell:")
    print_moments(moments)

    tracker += time.time()
    print("Runtime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
    exit(0)