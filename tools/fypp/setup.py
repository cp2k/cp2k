# -*- coding: utf-8 -*-
from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='fypp',

    version='1.2',

    description='Python powered Fortran preprocessor',
    long_description=long_description,

    url='https://bitbucket.org/aradi/fypp',

    author='Bálint Aradi',
    author_email='baradi09@gmail.com',

    license='BSD',

    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Code Generators',
        'Topic :: Software Development :: Pre-processors',

        'License :: OSI Approved :: BSD License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    keywords='fortran metaprogramming pre-processor',

    package_dir={'': 'src'},
    py_modules=['fypp'],

    entry_points={
        'console_scripts': [
            'fypp=fypp:run_fypp',
        ],
    },
)
