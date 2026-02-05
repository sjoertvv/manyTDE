# from distutils.core import setup
# setup(name='manyTDE',
#       version='1.0',
#       include_package_data = True
#       )

from setuptools import setup
import os.path
import json

exts =[] 

def read(rel_path):
      here = os.path.abspath(os.path.dirname(__file__))
      with open(os.path.join(here, rel_path), 'r', encoding='utf8') as fp:
            return fp.read()

def get_version():
      return json.load(open("manyTDE/data/sources/AT2019dsg.json",'rb'))["catalog_version"]
      raise RuntimeError("Unable to find version string.")    
  
    
if __name__ == "__main__" :
  setup(name='manyTDE',
        version=get_version(),
        author='Sjoert van Velzen',
        author_email='sjoert@strw.leidenuniv.nl',
        packages=['manyTDE', ],
        description='Optical TDEs with host information and light curves.',
        long_description=read('README.md'),
        install_requires=["emcee",
                          "numpy",
                          "astropy",
                          "matplotlib", 
                          "corner", 
                          ],
        ext_modules = exts,
        python_requires='>=3.5, <4'
      )
