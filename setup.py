from setuptools import setup, find_packages
import re

version = None
for line in open("./ArraySplitter/__init__.py"):
    m = re.search("__version__\s*=\s*(.*)", line)
    if m:
        version = m.group(1).strip()[1:-1]  # quotes
        break
assert version

setup(
    name='ArraySplitter',
    version=version,
    packages=["ArraySplitter"],
    package_data={"": ["README.md"]},
    python_requires='>=3.6',
    include_package_data=True,
    scripts=[],
    license="BSD",
    url='https://github.com/aglabx/ArraySplitter',
    author='Aleksey Komissarov',
    author_email='ad3002@gmail.com',
    description='De Novo Decomposition of Satellite DNA Arrays into Monomers within Telomere-to-Telomere Assemblies',
    install_requires=[
        'PyExp',
        'intervaltree',
        'editdistance',
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
