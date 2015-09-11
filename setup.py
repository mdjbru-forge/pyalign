# Written usingn resources from:
# https://packaging.python.org/en/latest/distributing.html#working-in-development-mode
# https://github.com/pypa/sampleproject/blob/master/setup.py

from setuptools import setup, find_packages

setup(name = "pyalign",
      version = "0.0.7",
      py_modules = ["pyalign",
                    "pyalignScripts"],
      install_requires = ["biopython",
                          "pygenes"],
      entry_points =  {
          "console_scripts" : [
              "pyalign=pyalignScripts:main"
          ]
      }
)
