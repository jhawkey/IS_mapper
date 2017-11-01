#!/usr/bin/env python

from distutils.core import setup

setup(
    name='ISMapper',
    version='0.1.5.1',
    author='Jane Hawkey',
    author_email='hawkey.jane@gmail.com',
    packages=['ismap'],
    scripts=['scripts/ismap.py', 'scripts/binary_table.py', 'scripts/compiled_table.py', 'scripts/create_genbank_table.py',
            'scripts/slurm_ismap.py', 'scripts/create_typing_out.py', 'scripts/slurm_ismap_sg.py'],
    entry_points={
        'console_scripts': ['ismap = ismap.ismap:main']
    },
    package_dir = {'ismap': 'scripts'},
    #package_data={'srst2': ['data/resistance.*']},
    url='http://jhawkey.github.io/IS_mapper/',
    license='LICENSE.txt',
    description='Discovering Insertion Sequences from Short Read Sequence Data',
    long_description=('This program is designed to take Illumina'
                      'sequence data, a query gene (such as an insertion sequence)'
                      'and a reference genome, either for finding novel positions in,'
                      'or an assembly where the query is either assembled or broken over'
                      'contigs.'),
    install_requires=[
        # Although we depend on scipy, which depends on numpy, we don't
        # specify the dependencies here because they don't play well with
        # any Python installing system, such as pip or easy_install.
        # So we assume the user has already installed the dependencies
        # themselves.
        #"numpy >= 1.7.1",
        #"scipy >= 0.12.0",
        #"biopython == 1.63",
    ],
)
