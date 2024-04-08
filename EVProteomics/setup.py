from setuptools import setup, find_packages

setup(
    name='EVProteomics',
    version='0.1',
    description='A collection of tools to help with proteomic analysis of EVs (extracellular vesicles)',
    author='Sydney DAmaddio',
    author_email='sydney.damaddio@wyss.harvard.edu',
    packages=find_packages(),
    install_requires=[
        'biolib',
        'gzip',
        'os',
        'requests',
        'pathlib',
        'scipy',
        'matplotlib',
        'matplotlib.pyplot',
        'numpy',
        'pandas',
        'requests',
        'seaborn',
        'io',
    ]
)