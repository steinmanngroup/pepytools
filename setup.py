import sys
#from distutils.core import setup
from numpy.distutils.core import Extension, setup

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


def setup_pepylib():

    setup(
          name="pepylib",
          version="0.1",
          author="Casper Steinmann",
          author_email="casper.steinmann@gmail.com",

          package_dir={'pepylib': 'src'},
          packages=['pepylib'],
          scripts=['bin/pepy_add'],
          ext_package = 'pepylib',
          ext_modules = [ext_field, ext_intersect, ext_qmfields],

)

if __name__ == '__main__':
    setup_pepylib()
