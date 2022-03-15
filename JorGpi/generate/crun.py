# -*- coding: utf-8 -*-

from os import scandir,mkdir
import re
from ctypes import CDLL,c_char_p
import importlib
import setuptools

class Crun:
    compilerOptions= { 'include_dirs'       : ['/usr/local/include'],
                       'extra_compile_args' : ['-O3','-Wall','-pedantic']}

    def __init__(self,*args,**kwargs):
        importlib.reload(setuptools)
        try:
            mkdir('build')
        except Exception as ext:
            print(ext,": ignoring")
        self.files = args
        self.compilerOptions.update(kwargs)

        if 'name' in kwargs:
            self.name = kwargs['name']
            del self.compilerOptions['name']
        else:
            self.name = re.sub('^.*/','',args[0])
            self.name = re.sub('\..*$','',self.name)

        self.module = setuptools.Extension(self.name, sources = list(args), **self.compilerOptions)

        setuptools.setup(
                script_args=['build','--build-lib=%s-lib'%(self.name)],
                ext_modules=[self.module])

        setuptools.setup(name=self.name,script_args=['clean'],
               ext_modules = [self.module])

        self.search_for_library()

    def search_for_library(self):
        match = [re.match('^%s.*\.so'%self.name,entry.name)\
                 for entry in scandir('%s-lib'%(self.name))]
        try:
            while None in match:
                match.remove(None)
        except ValueError as err:
            if "x not in list" not in str(err):
                print('ValueError: ',str(err))
                exit(-1)
        try:
            self.library = CDLL('%s-lib/%s'%(self.name,match[0].string))
        except UnboundLocalError:
            print("I couldn't find directory")
            exit(-1)

    def __call__(self,function,*args):
        argStr = []
        for arg in args:
            if isinstance(arg,str):
                argStr.append(c_char_p(arg.encode('utf-8')))
            else:
                argStr.append(arg)
        f = eval("self.library.%s"%function)
        return f(*argStr)

    def __del__(self):
        setuptools.setup(script_args=['clean','--all','--build-lib=%s-lib'%(self.name)],
               ext_modules = [self.module])
        importlib.reload(setuptools)
