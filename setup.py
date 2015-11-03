import sys
#from distutils.core import setup
from numpy.distutils.core import Extension, setup

from src import __version__, __license__, __author__, __copyright__, __email__, __description__, __url__

# setup fortran extensions
ext_field = Extension(name = 'field',
                      sources = ['src/field.f90'],
                      language = 'fortran',
                      extra_f90_compile_args = ['-fopenmp'],
                      extra_link_args = ['-lgomp'])

ext_intersect = Extension(name = 'intersect',
                          sources = ['src/intersect.f90'],
                          language = 'fortran')

ext_qmfields = Extension(name = 'qm_fields',
                         sources = ['src/qm_fields.f90'],
                         language = 'fortran')

# use README.md as long description
def readme():
    with open('README.md') as f:
        return f.read()

def setup_pepylib():

    setup(
          name="pepylib",

          # metadata
          version=__version__,
          author=__author__,
          author_email=__email__,
          license = __license__,
          platforms = 'Any',
          description = __description__,
          long_description = readme(),
          keywords = 'polarizable embedding potentials',
          url = __url__,

          # set up package contents
          package_dir={'pepylib': 'src'},
          packages=['pepylib'],
          scripts=['bin/pepy_add'],
          ext_package = 'pepylib',
          ext_modules = [ext_field, ext_intersect, ext_qmfields],
)

if __name__ == '__main__':
    setup_pepylib()
