import setuptools

#with open('README.md', 'r') as fh:
#   long_description = fh.read()

setuptools.setup(
    name='wg-blimp',
    version='0.9.6',
    author='Marius Woeste',
    author_email='mar.w@wwu.de',
    description='WGBS methylation analysis pipeline',
    packages=[
        'cli',
        'snakemake_wrapper',
        'shiny'
    ],
    py_modules=[
        'cli'
    ],
    package_data={
        'snakemake_wrapper': [
            'annotation/*',
            'scripts/*.R',
            'submodules/camel/**/*',
            'Snakefile'
        ],
        'shiny': [
            '*.R'
        ],
        'cli': [
            'annotation_download_links.yaml',
            'optionals.yaml'
        ]
    },
    install_requires=[
        'click',
        'ruamel.yaml',
        'snakemake'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research'
    ],
    entry_points='''
        [console_scripts]
        wg-blimp=cli.cli:main
    '''
)
