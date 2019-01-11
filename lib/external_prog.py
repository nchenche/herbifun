# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 18:17:17 2018

@author: nicolas.chevrollier
"""

import os

def run_hmmbuild(clwfile=None, hmmname=None, output=None, logger=None):
    if clwfile is not None:
        if output is not None:
            if hmmname is None:
                cmd = 'hmmbuild '+'-o /dev/null '+output+' '+clwfile
                logger.info('Running '+ cmd)
                os.system(cmd)
            else:
                cmd = 'hmmbuild '+'-o /dev/null -n '+hmmname+' '+output+' '+clwfile
                logger.info('Running '+ cmd)
                os.system(cmd)
        else:
            print('HMM output name is not defined'.format())
    else:
        print('HMM input name is not defined'.format())


def run_hmmsearch(hmmfile=None, seqdb=None, output=None, logger=None):
    if hmmfile is not None:
        if output is not None:
            if seqdb is not None:
                cmd = 'hmmsearch'+' -o /dev/null '+'--domtblout '+output+' '+hmmfile+' '+seqdb
                logger.info('Running '+ cmd)
                os.system(cmd)
            else:
                print('No sequences target provided...'.format())
        else:
            print('No output name defined...'.format())
    else:
        print('No hmmfile provided...'.format())


def run_usearchclust(msafile=None, output=None, identity=0.9, logger=None):
    """
    - Usearch clustering to remove sequence redundancy
    - Writes the non-redundant output file as _nr.mfs
    """
    if output is None:
        output = msafile.split('.')[0]+'_nr.msa'
    else:
        cmd = 'usearch -sort length -cluster_fast '+msafile+' -id '+str(identity)+' -centroids '+output+' -quiet'
        logger.info('\nRunning {}'.format(cmd))
        os.system(cmd)


def run_muscle(msafile=None, output=None, logger=None):
    """
    - Multiple alignment with muscle
    """
    if msafile is not None:
        if output is None:
            output = msafile.split('.')[0]+'.clw'
        else:
            cmd = 'muscle -in '+msafile+' -out '+output+' -clw -quiet'
            logger.info('Running '+ cmd)
            os.system(cmd)
    else:
        print('A fasta file is required...'.format())


def hmmscan(proteome_filename=None, scanoutput=None, hmmdb_filename=None):
    """
    - Calls hmmscan from hmmer
    - Returns outputfile in domtblout format
    """
    cmd = 'hmmscan '+'-o /dev/null '+'--domtblout '+scanoutput+' '+hmmdb_filename+' '+proteome_filename

    print(cmd)
    os.system(cmd)
