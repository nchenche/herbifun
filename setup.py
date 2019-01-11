from distutils.core import setup
#from setuptools import setup
setup(name='HMMbuilder',
	packages = ['lib'],
	scripts = ['bin/hmmbuilder.py', 'bin/annotater.py', 'bin/usearch'],

      	data_files=[('datas', ['datas/A.msa', 'datas/AT.msa', 'datas/KS.msa', 'datas/PP.msa',
		'datas/database.hmm', 'datas/database.hmm.h3f', 'datas/database.hmm.h3i', 'datas/database.hmm.h3m', 'datas/database.hmm.h3p',
		'datas/annotation.rules',
		'datas/mgg_70-15_8.fasta','datas/mgg_70-15_8.domtblout',
		'requirements.txt'])])


import os
import sys

if sys.argv[1] != 'sdist':
    cmd = 'pip install -r requirements.txt'
    print('\n{}'.format('-'*25))
    print('{}'.format(cmd))
    os.system(cmd)
    print('\n{}\n'.format('-'*25))

