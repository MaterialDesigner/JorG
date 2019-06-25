#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
from POSCARloader import POSCARloader
from loadsave import save_vanilla_POSCAR

loader=POSCARloader('POSCAR')
loader.parse()
a = loader()
print(a)
save_vanilla_POSCAR('POSCAR_fixed',a)
