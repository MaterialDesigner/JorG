# -*- coding: utf-8 -*-
from sys import path
path.insert(0,r'../')

import numpy as np
from aux.periodic import maskFull

def find_unique_flips(referenceAtom,crystal,symmetry,cutOff,
                      atomTypeMask=maskFull,
                      Wyckoffs="abcdefghijklmnopqrstuvwxyz",
                      wyckoffDict={'a':'a','b':'b','c':'c','d':'d','e':'e','f':'f','g':'g','h':'h','i':'i','j':'j','k':'k','l':'l','m':'m','n':'n','o':'o','p':'p','q':'q','r':'r','s':'s','t':'t','u':'u','v':'v','w':'w','x':'x','y':'y','z':'z'}, 
                      logAccuracy=2):

    distances = []
    for i,atom in enumerate(crystal):
       d = np.around(np.linalg.norm(atom[1]-referenceAtom[1]), decimals=logAccuracy)
       if ((d,atom[0]) not in distances and
            wyckoffDict[symmetry['wyckoffs'][i]] in Wyckoffs and
            d > 0.0 and d <= cutOff and '$'+atom[0]+'$' in atomTypeMask):
            distances.append((d,atom[0]))
    
    distances = np.array(distances, dtype=[('distance', np.float), ('element', 'U3')])
    distances.sort(order='distance')
    
    checker = [ False        for x in range(len(crystal))] 
    flipper = [] 
    case = 1
    for (d,n) in distances:
        for i,atom in enumerate(crystal):
            if checker[i]:
                continue
            if d == np.around(np.linalg.norm(atom[1]-referenceAtom[1]),decimals=logAccuracy)\
              and atom[0] == n:
                for j in np.argwhere(symmetry['equivalent_atoms'] == symmetry['equivalent_atoms'][i]).flatten():
                    if d == np.around(np.linalg.norm(crystal[j][1]-referenceAtom[1]),decimals=logAccuracy)\
                      and wyckoffDict[symmetry['wyckoffs'][j]] in Wyckoffs:
                        checker[j] = True
                flipper.append((i,atom,d,wyckoffDict[symmetry['wyckoffs'][i]]))
                case += 1
    flipper = np.array(flipper, dtype=[('id', int), ('atom', list), ('distance', np.float), ('wyckoff', 'U1')])
    flipper.sort(order='distance')
    return flipper

#
#
#
#
#

def find_all_flips(referenceAtom,crystal,symmetry,cutOff,
                   atomTypeMask=maskFull,
                   Wyckoffs="abcdefghijklmnopqrstuvwxyz",
                   wyckoffDict={'a':'a','b':'b','c':'c','d':'d','e':'e','f':'f','g':'g','h':'h','i':'i','j':'j','k':'k','l':'l','m':'m','n':'n','o':'o','p':'p','q':'q','r':'r','s':'s','t':'t','u':'u','v':'v','w':'w','x':'x','y':'y','z':'z'}, 
                   logAccuracy=2):
   
    flipper = []
    for i,atom in enumerate(crystal): 
        if "$"+atom[0]+"$" in atomTypeMask                         \
            and wyckoffDict[symmetry['wyckoffs'][i]] in Wyckoffs   \
            and np.linalg.norm(atom[1]-referenceAtom[1]) > 1e-3:
            flipper.append(i)

#            and np.linalg.norm(atom[1]-referenceAtom[1]) <= cutOff 

    return flipper


def find_all_distances(reference, 
                       crystal8,
                       cutOff, flipper,
                       logAccuracy = 2):
    acceptableFlips = np.append(flipper, reference)

    size = int(len(crystal8)/8)
    distances = []
    for i in range(8):
        for flipA in acceptableFlips:
            for j in range(8):
                for flipB in acceptableFlips:
                    d = np.round(np.linalg.norm(crystal8[flipA+size*i][1] - crystal8[flipB+size*j][1]),logAccuracy)
                    if d not in distances \
                       and d <= cutOff:
                        distances.append(d)
    distances = np.array(distances)
    distances.sort()
    return distances
            