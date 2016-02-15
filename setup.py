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

def setup_pepytools():

    setup(
          name="pepytools",

          # metadata
          version=__version__,
          author=__author__,
          author_email=__email__,
          license = __license__,
          platforms = 'Any',
          description = __description__,
          long_description = readme(),
          keywords = 'polarizable embedding potential',
          url = __url__,

          # set up package contents
          package_dir={'pepytools': 'src'},
          packages=['pepytools'],
          scripts=['bin/pepy_add', 'bin/pepy_to_csv'],
          ext_package = 'pepytools',
          ext_modules = [ext_field, ext_intersect, ext_qmfields],
)

if __name__ == '__main__':
    setup_pepytools()
