from distutils.core import setup, Extension

module1 = Extension("nwchem", 
                    libraries = ['nwapi'],
                    library_dirs = ['./'],
                    sources = ["nwchem_wrap.c"])
setup(name="nwchem", version="1.0",
      ext_modules=[module1])
