
from setuptools import find_packages, setup

setup(name='conduction',
      version='1.0',
      description='A simulation package to measure the thermal transport properties of a carbon nanotube - polymer '
                  'composite using random walks. Developed in the Mullen theoretical condensed matter physics '
                  'research group at the University of Oklahoma (OU)',
      author='Timothy A. Burt',
      author_email='taburt@ou.edu',
      url='https://github.com/tab10/conduction',
      packages=find_packages(),
      install_requires=['numpy', 'matplotlib', 'scipy', 'mpi4py'], )
