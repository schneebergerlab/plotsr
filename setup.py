from setuptools import setup, Extension
from plotsr import __version__

setup(
    name="plotsr",
    version='{}'.format(__version__),
    packages=["plotsr",],
    license='MIT License',
    long_description=open('README.md').read(),
)