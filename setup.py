from setuptools import setup

with open("README.rst", 'r') as f:
    long_description = f.read()

setup(
   name='circosmetaplots',
   version='0.1',
   description='Visualizes metagenomics data for every other week over a year, with weather data',
   license="MIT",
   long_description=long_description,
   author='mjlundkvist',
   author_email='molu0019@ad.umu.se',
   url="http://github.com/mjlundkvist/circosmetaplots",
   packages=['circosmetaplots'],  #same as name
   install_requires=['pandas', 'numpy','subprocess', 'shutil'], #external packages as dependencies
  
)