from setuptools import setup, find_packages

setup(
  name = 'QAlign',
  version = '0.0.1',
  description = 'QAlign: Aligning nanopore reads accurately using current-level modeling',
  url='https://github.com/olivomao/QAlign',
  author='Dhaivat Joshi and Shunfu Mao',
  install_requires=[
    'numpy',
  ],
  packages = find_packages()
)