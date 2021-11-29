import re
from setuptools import setup
from os import makedirs,environ
from subprocess import call
import shlex
from sys import argv,path,version_info

packages_jorgpi = ['JorGpi',
                   'JorGpi.generate',
                   'JorGpi.pickup',
                   'JorGpi.utilities',
                   'JorGpi.utilities.KPOINTS',
                   'JorGpi.utilities.fixPOSCAR',
                   'JorGpi.aux',
                   'JorGpi.geometry']
executables_jorgpi = ['JorGpi/bin/JorGpi-demagnetize',
                      'JorGpi/bin/JorGpi-KPOINTS',
                      'JorGpi/bin/JorGpi-pickup',
                      'JorGpi/bin/JorGpi-POSCAR',
                      'JorGpi/bin/JorGpi-startup']

requirements_jorgpi = []

VERSION='0.2.1'#.heisenberg'

if __name__ == '__main__':
    try:
        myself=environ["PWD"]
    except KeyError:
        print("PWD not set! Try running sudo -E python3 setup.py install instead.")
    options=""
    try:
        if argv[1] != 'install':
            exit(-1)
    except IndexError:
        exit(-1)
    for arg in argv[2:]:
        options+=" %s"%arg

    options=re.sub('=',' ',options)
    options=re.sub('~',environ['HOME'],options)
    if '--install_reqs' in argv:
        print('Installing requirements')
        call(shlex.split('./install-requirements.sh'+options),shell=False)
        argv.remove('--install_reqs')

    if '--user' in argv:
        try:
            locDir = environ['PYTHONUSERBASE']
            with open(environ['HOME']+'/.bashrc','a+') as bashrc:
                bashrc.write('\nexport PYTHONPATH=%s/lib/python%d.%d/site-packages/JorGpi-%s-py%d.%d.egg:$PYTHONPATH\n'%(locDir,version_info[0],version_info[1],VERSION,version_info[0],version_info[1]))
            makedirs(locDir,exist_ok=True)
        except KeyError:
            locDir = environ['HOME']+'/.local'
            print('Installing in %s'%(environ['HOME']))
        options+=' --prefix %s'%locDir
    try:
        with open('requirements.txt','r') as reqs:
            for line in reqs.readlines():
                package = re.sub('[><=]+.*','',line)
                package = re.sub('\s','',package)
                if len(package)>0:
                    requirements_jorgpi.append(package)
    except FileNotFoundError:
        requirements_jorgpi = [ 'numpy', 'scipy', 'spglib', 'matplotlib', 'setuptools', 'defusedxml' ]
    setup(name='JorGpi',
          version=VERSION,
          description='JorGpi DFT-to-Heisenberg mapping module',
          long_description="""
          JorGpi DFT-to-Heisenberg mapping module
          """,
          author='Andrzej P. Kądzielawa',
          author_email='andrzej.piotr.kadzielawa@vsb.cz',
          url='https://github.com/Mellechowicz/JorG',
          packages=packages_jorgpi,
          install_requires=requirements_jorgpi,
          provides=['JorGpi'],
          scripts=executables_jorgpi,
          platforms=['POSIX'],
          license='',
          test_suite='nose.collector',
          tests_require=['nose'])

    call(shlex.split('./install-gsl.sh'+options),shell=False)
    call(shlex.split('./install-asa.sh'+options),shell=False)
