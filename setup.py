#!/usr/bin/env python3

import setuptools

setuptools.setup(
    name="ISMapper",
    version="2.0",
    author="Jane Hawkey",
    description="Detection of insertion sequence locations from short read sequence data",
    url="https://github.com/jhawkey/IS_mapper",
    packages=['ismap'],
    scripts=['scripts/ismap.py', 'scripts/compiled_table.py', 'scripts/create_output.py', 'scripts/mapping_to_query.py',
             'scripts/mapping_to_ref.py', 'scripts/read_grouping.py', 'scripts/run_commands.py'],
    entry_points={
        'console_scripts': ['ismap = ismap:main']
    },
    package_dir = {'ismap': 'scripts'},
    license='LICENSE.txt',
)