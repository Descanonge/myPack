
from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='mypack-descanonges',
      version='1.0',
      description='Convenience functions',
      long_description=readme(),
      url='http://github.com/Descanonges/myPack',
      author='Clément HAËCK',
      author_email='clement.haeck@postes.net',
      license='',
      packages=['myPack'],
      install_requires=['regulargrid'])
