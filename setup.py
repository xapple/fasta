from setuptools import setup, find_packages

setup(
    name             = 'fasta',
    version          = '2.0.4',
    description      = 'The fasta python package enables you to deal with '
                       'biological sequence files easily.',
    license          = 'MIT',
    url              = 'http://github.com/xapple/fasta/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = find_packages(),
    install_requires = ['plumbing>=2.8.1', 'autopaths>=1.4.2',
                        'biopython', 'numpy', 'sh', 'tqdm'],
    long_description = open('README.md').read(),
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)
