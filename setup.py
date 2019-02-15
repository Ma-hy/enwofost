import os
from setuptools import setup

try:
    version = os.environ['SIAC_VERSION']
except:
    version = '0.0.3'

with open('README.md', 'rb') as f:
    readme = f.read().decode()

setup(name                          = 'enwofost',
      version                       = version,
      description                   = 'Ensembles of wofost',
      long_description              = readme,
      long_description_content_type ='text/markdown',
      author                       = 'Hongyuan Ma',
      author_email                 = 'Hongyuan.Ma@ucl.ac.uk',
      classifiers                  = ['Development Status :: 4 - Beta',
                                      'Programming Language :: Python :: 2.7',
                                      'Programming Language :: Python :: 3.6'],
      install_requires             = ['gdal>=2.1', 'numpy>=1.13', 'scipy>=1.0', 'pcse',
                                      'netCDF4', 'cdsapi'],
      url                          = 'https://github.com/Ma-hy/enwofost',
      license                      = "GNU Affero General Public License v3.0",
      include_package_data         = True,
      packages                     = ['enwofost'],
     )
