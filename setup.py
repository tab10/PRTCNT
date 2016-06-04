
from setuptools import find_packages, setup

setup(name='conduction',
      version='0.1',
      description='Package to study nanotubes in a medium',
      author='Timothy Burt',
      author_email='timcantango@gmail.com ',
      url='https://github.com/tab10/conduction',
      packages=find_packages(),
      install_requires=['numpy','matplotlib'],)