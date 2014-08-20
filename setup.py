from distutils.core import setup

setup(
        name             = 'fasta',
        version          = '1.0.0',
        description      = 'The fasta python package enables you to deal with biological sequence files easily',
        long_description = open('README.txt').read(),
        license          = 'MIT',
        url              = 'http://xapple.github.com/fasta/',
        author           = 'Lucas Sinclair',
        author_email     = 'lucas.sinclair@me.com',
        classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
        packages         = ['fasta'],
        install_requires = ['sh', 'biopython'],
    )
