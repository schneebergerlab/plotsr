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
    packages=["plotsr"],
    install_requires=[
        'matplotlib  ==3.3.4',
        'numpy ==1.21.2',
        'pandas ==1.2.4',
    ],
    py_modules=["plotsr.func", "plotsr.main"],
    scripts=['bin/plotsr'],
    long_description=open('README.rst').read(),
)
