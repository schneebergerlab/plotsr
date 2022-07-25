from setuptools import setup, Extension
from plotsr import __version__

setup(
    name="plotsr",
    version='{}'.format(__version__),
    description='Package to plot structural annotations between multiple genomes',
    author='Manish Goel',
    author_email='mnshgl0110@gmail.com',
    url='https://github.com/schneebergerlab/plotsr',
    license='MIT License',
    license_files=('LICENSE',),
    packages=["plotsr", "plotsr.scripts"],
    py_modules=["plotsr.scripts.plotsr",
                "plotsr.scripts.func"],
    scripts=['bin/plotsr'],
    long_description=open('README.rst').read(),
)
