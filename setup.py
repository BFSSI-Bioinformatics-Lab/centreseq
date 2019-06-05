import setuptools

__version__ = '0.2.0'
__author__ = ['Forest Dussault', 'Adrian Verster', 'Nicholas Petronella']
__email__ = 'forest.dussault@canada.ca'

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="centreseq",
    install_requires=['click>=7.0',
                      'pandas>=0.23.0',
                      'dataclasses>=0.6',
                      'xlsxwriter>=1.1.8',
                      'tqdm>=2.2.3',
                      'biopython>=1.70',
                      'scipy>=1.1',
                      'pycurl>=7.43',
                      'cyvcf2>=0.11.3'
                      ],
    python_requires='~=3.6',
    description="Fast generation of a core genome among bacterial strains",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bfssi-forest-dussault/centreseq",
    include_package_data=True,  # Must be supplemented by MANIFEST.in file containing paths to extra files
    install_package_data=True,
    packages=setuptools.find_packages(),
    version=__version__,
    author=__author__,
    author_email=__email__,
    entry_points={
        'console_scripts': [
            'centreseq=centreseq.centreseq:cli'
        ]
    }
)
