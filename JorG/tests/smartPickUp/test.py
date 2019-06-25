#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import argv,path
path.insert(0,r'../../')
import time
import numpy as np
from POSCARloader import POSCARloader
from pickup.pickup import SmartPickUp
from itertools import product

def main(**args):
    pass

if __name__ == '__main__':
    tracker  = -(time.time())

    POSCARs = [ "%s/POSCAR"%arg for arg in argv[1:] ]
    if not POSCARs:
        print("No files provided")
        argv.append("../../_VASP/Fe/noFlip")
        argv.append("../../_VASP/Fe/flip00000")
        POSCARs = [ "%s/POSCAR"%arg for arg in argv[1:] ]

    print(*POSCARs)
    loader   = POSCARloader(*POSCARs,spam=False)
    loader.parse()

    print("Running for NN=1, \'Fe\':")
    pickerUpper = SmartPickUp(1,'Fe')
    pickerUpper.read(*argv[1:])
    for unit in ['eV', 'meV', 'Ry', 'mRy', 'He', 'mHe', 'K']:
        print("Exchange interaction in %s:"%unit)
        Js = pickerUpper.solve(units=unit)
        for i,typeName in enumerate(pickerUpper.types):
            print(("  %s:\t"+len(Js[0])*"% 11.7f ")%(typeName,*Js[i],))

    tracker += time.time()
    print("Runntime of %02d:%02d:%02d.%09d"%(int(tracker/3600),int(tracker/60),int(tracker),int(1e9*tracker)))
