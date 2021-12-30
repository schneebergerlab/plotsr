from setuptools import setup, Extension
from plotsr import __version__

setup(
    name="plotsr",
    version='{}'.format(__version__),
    packages=["plotsr"],
    py_modules=["plotsr.func"],
    license='MIT License',
    scripts=['plotsr/plotsr',],
    long_description=open('README.md').read(),
)
