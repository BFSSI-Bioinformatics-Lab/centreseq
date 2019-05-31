import setuptools

__version__ = '0.1.0'
__author__ = ['Forest Dussault', 'Adrian Verster', 'Nicholas Petronella']
__email__ = '0.1.0'

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="centreseq",
    install_requires=['click', 'pandas', 'dataclasses', 'xlsxwriter', 'tqdm', 'biopython'],
    python_requires='~=3.6',
    description="Package for checking for the presence/absence of markers against a set of samples",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bfssi-forest-dussault/centreseq",
    package_data={'centreseq': ['*']},
    install_pacakage_data=True,
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
