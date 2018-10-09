# -*- coding: utf-8 -*-
import argparse as ap
import numpy as np
import aux.periodic as periodic

class options:
    keys = ["cutOff", "neighbor", "Wyckoffs", "reference", "input", "output", "mask", "symmetry"]
    def __init__(self, *args):
        self.parser = ap.ArgumentParser(description='Find minimal number of unique spin-flips')
        typeOfRange = self.parser.add_mutually_exclusive_group()
        typeOfRange.add_argument('--cutOff', '-R', default=-1.0, type=np.float,
                                 help='a cut-off distance (in Å) for calculations')
        typeOfRange.add_argument('--neighbor', '-N', default=-1, type=int,
                                 help='a rank of the last Neighbor taken into account')
        self.parser.add_argument('--Wyckoffs', '-W', default='abcdefghijklmnopqrstuvwxyz',
                                 help='narrows down the atomic selection to the atoms in positions defined by string (eg. \'abc\')')
        self.parser.add_argument('--symmetry', '-S', action='store_true',
                                 help='symmetry run only (default False)')
        self.parser.add_argument('--spin-orbit', '--SOC', action='store_true',
                                 help='is sping-orbit coupling enabled (default False)')
        self.parser.add_argument('--reference', '-r', default=1, type=int,
                                 help='number of reference atom in inputFile')
        self.parser.add_argument('--input', '-i', default='POSCAR',
                                 help='input POSCAR file')
        self.parser.add_argument('--output', '-o', default='output/POSCAR',
                                 help='output POSCAR file')
        self.parser.add_argument('--elements','-E', 
                                 help='string of all elements taken into account (eg. \'CuO\')')
        self.parser.add_argument('--group',choices=range(1,19), type=int, nargs='+',
                                 help='group number (eg. 1 <=> \'HLiNaKRbCsFr\')')
        self.parser.add_argument('--period',choices=periodic.periodNames, nargs='+',
                                 help='period name (eg. 2d <=> \'%s\')'%periodic.mask3d)
        self.parser.add_argument('--block',choices=periodic.blockNames, nargs='+',
                                 help='block name (eg. P <=> \'%s\')'%periodic.maskP)

        self.opt = self.parser.parse_args(args[1:])

    def __str__(self):   
        from textwrap import wrap
        output = ""
        for opt in self.opt.__dict__:
            output += "\n".join(wrap("%s =  %s"%("{:<10}".format(str(opt)),str(self.opt.__dict__[opt])),70,subsequent_indent=14*" "))
            output += "\n"
        return output[:-1]
    
    def __call__(self, key):
        if key not in self.keys:
            print("No key %s defined, please try: "%key)
            output = ""
            for k in self.keys:
                output += k + ", "
            print(output[:-2])
            exit(-300)
        elif key == 'reference':
            return self.opt.__dict__['reference'] - 1; 
        elif key == 'neighbor':
            neighbor = self.opt.__dict__['neighbor'] 
            if self.opt.__dict__['neighbor'] < 0:
                if self.opt.__dict__['cutOff'] < 0.0:
                    return 1
                else:
                    print("No neighbor defined: use cutOff")
                    return None                    
            return neighbor 
        elif key == 'cutOff':
            cutOff = self.opt.__dict__['cutOff'] 
            if cutOff < 0.0:
                print("No cutOff defined: use neighbor")
                return None                    
            return cutOff    
        elif key in self.opt.__dict__:
            return self.opt.__dict__[key]
        elif key == 'mask':
            output = ''
            if self.opt.__dict__['period'] is not None:
                for el in self.opt.__dict__['period']:
                    output += eval("periodic.mask"+el)
            if self.opt.__dict__['group'] is not None:
                for el in self.opt.__dict__['group']:
                    output += eval("periodic.maskG%d"%el)
            if self.opt.__dict__['block'] is not None:
                for el in self.opt.__dict__['block']:
                    output += eval("periodic.mask"+el)
            if self.opt.__dict__['elements'] is not None: 
                output += self.opt.__dict__['elements']
            if output == '':
                return periodic.maskFull
            return output
           
                
