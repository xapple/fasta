from distutils.core import setup

setup(
    name             = 'fasta',
    version          = '1.0.2',
    description      = 'The fasta python package enables you to deal with biological sequence files easily.',
    license          = 'MIT',
    url              = 'http://github.com/xapple/fasta/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = ['fasta'],
    install_requires = ['plumbing', 'sh', 'biopython'],
    long_description = open('README.md').read(),
)
