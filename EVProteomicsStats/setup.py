from setuptools import setup, find_packages

setup(
    name='mypackage',
    version='0.1',
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
        'io'
    ],
    entry_points={
        'console_scripts': [
            'mypackage_script = mypackage.module1:main',
        ],
    },
)