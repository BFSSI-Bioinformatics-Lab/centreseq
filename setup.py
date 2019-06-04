import setuptools

__version__ = '0.1.3'
__author__ = ['Forest Dussault', 'Adrian Verster', 'Nicholas Petronella']
__email__ = 'forest.dussault@canada.ca'

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="centreseq",
    install_requires=['click', 'pandas', 'dataclasses', 'xlsxwriter', 'tqdm', 'biopython'],
    python_requires='~=3.6',
    description="Core genome finder",
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
