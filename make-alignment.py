from distutils.core import setup, Extension

alignment = Extension('alignment',
                    sources = ['alignment.cpp'],
                    define_macros=[('__PYTHON__', None)])

setup (name = 'alignment',
       version = '0.1',
       description = "Wrapper for Conan's alignment.cpp code",
       ext_modules = [alignment])
