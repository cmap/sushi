# docker/prism-tools/setup.py
from setuptools import setup, find_packages

setup(
  name="prism-tools",        # distribution name, can have dashes
  version="0.1.0",
  packages=find_packages(),  # will discover prism_tools
)
