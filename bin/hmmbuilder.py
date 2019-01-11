#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 15:45:57 2018

@author: nicolas.chevrollier
"""
#==============================================================================
#  Imports
#==============================================================================
import os
import sys
import shutil
import datetime
date = datetime.datetime.now()
date = '_'.join(str(date).split('.')[0].split()).replace(':','-')

#import _path2lib
#_path2lib.link_bin2lib()

import lib.external_prog as external_prog
import lib.seqio as seqio
import lib.hmmdomtbl as hmmdomtbl
import lib.logger as logger

import argparse

#==============================================================================
#  Functions
#==============================================================================
def buildhmm(logger=None,input_msa=None, proteome_filename=None, path2output=None, iter_nb=1, clust=True,
         identity = 0.9, coverage_co=0.0, cval_co=0.01, ival_co=0.01, acc_co=0.6):  
    
    txt = '{}'.format(40*'-')
    txt = txt + ' Running hmmbuilder iteration {} '.format(iter_nb)
    if iter_nb == 1:
        txt = txt + '{}'.format(40*'-')
    else:
        txt = txt + '{}\n'.format(40*'-')
    logger.info(txt)
    # Inputs
    input_msaraw = input_msa
    
    # gets pattern name of output files
    sepdir = os.path.dirname(input_msa)
    if len(sepdir) != 0:
        dom_name = input_msa.split(sepdir)[-1].split('/')[-1].split('.')[0]
    else:
        dom_name = input_msa.split('.')[0]
    
    
    # output path    
    if not os.path.exists(path2output):
        os.makedirs(path2output)
    
    
    # output filenames
    output_msanr = path2output+dom_name+'-'+str(iter_nb)+'_nr.msa'
    output_clw = path2output+dom_name+'-'+str(iter_nb)+'_nr.clw'
    output_hmm = path2output+dom_name+'-'+str(iter_nb)+'_nr.hmm'
    output_domtbl = path2output+dom_name+'-'+str(iter_nb)+'_nr.domtblout'
    output_msanew = path2output+dom_name+'-'+str(iter_nb+1)+'_new.msa'
    output_msahybrid = path2output+dom_name+'-'+str(iter_nb+1)+'_hybrid.msa'
    output_msahybridnr = path2output+dom_name+'-'+str(iter_nb+1)+'_nr.msa'
    
    
    if clust == True:
        # Removes redundancy
        external_prog.run_usearchclust(msafile=input_msaraw, output=output_msanr, identity=identity, logger=logger)
    
    # Aligns sequences in msa file
    external_prog.run_muscle(msafile=output_msanr, output=output_clw, logger=logger)
    
    # hmmbuild
    external_prog.run_hmmbuild(clwfile=output_clw, hmmname=dom_name, output=output_hmm, logger=logger)
    
    # hmmsearch
    external_prog.run_hmmsearch(hmmfile=output_hmm, seqdb=proteome_filename, output=output_domtbl, logger=logger)
    
    
    # parsing hmmsearch output
    txt = 'Parsing {} searching for domains meeting the following criteria:\n'.format(output_domtbl)
    txt = txt + '  - envlen >= qlen * {} & dom_cval <= {} & dom_ival <= {} & acc >= {}\n'.format(coverage_co, cval_co, ival_co, acc_co)   
    logger.info(txt)
    hmmsearch_hits = hmmdomtbl.domtbl_dic(domtblout=output_domtbl, domtype='search')
    hmmdomtbl.filtering_hits(domtbl=hmmsearch_hits, coverage_co=coverage_co, cval_co=cval_co, ival_co=ival_co, acc_co=acc_co, logger=logger)
    
    # gets sequences
    proteome_fasta = seqio.read_fasta(sequences=proteome_filename)  
    hmmdomtbl.get_domtblseq(domtbl_dic=hmmsearch_hits, fasta=proteome_fasta)
    
    # writes hits fasta
    dom_fasta = hmmdomtbl.get_domtblfasta(hmmsearch_hits)
    seqio.write_fasta(seqfasta=dom_fasta, output=output_msanew)
    
    # concatenates new hits fasta with the original input fasta    
    seqio.concat_file(output_msahybrid, [output_msanew, input_msa])
    
    # Removes redundancy
    external_prog.run_usearchclust(msafile=output_msahybrid, output=output_msahybridnr, identity=identity, logger=logger)
    
   
    return output_msahybridnr, output_clw, output_hmm, hmmsearch_hits


#==============================================================================
#  Argument parsing
#==============================================================================
parser = argparse.ArgumentParser(description='Iterative building of hmm profiles')
parser.add_argument("-seqdb", required = True, nargs = "?", help = "Sequences used to learn hmm profile (fasta format)")
parser.add_argument("-fasta", required = True, nargs = "?", help = "Sequence(s) used as first seed (fasta format)")
parser.add_argument("-dir", required = True, nargs = "?", help = "Output directory")
parser.add_argument("-identity", required = False, help = "Sequence identity threshold to remove redundancy in seeds'sequences", type=float)
parser.add_argument("-cov", required = False, help = "Minimum percentage of coverage alignment between hmm hit and hmm profile (0.0)", type=float)
parser.add_argument("-cval", required = False, help = "hmmer conditional e-value cutoff (0.01)", type=float)
parser.add_argument("-ival", required = False, help = "hmmer independant e-value cutoff (0.01)", type=float)
parser.add_argument("-acc", required = False, help = "hmmer mean probability of the alignment accuracy between each residues of the target and the corresponding hmm state (0.6)", type=float)
args = parser.parse_args()

# Optional parameters
if args.identity is not None:
    identity = args.identity
else:
    identity = 0.9
if args.cov is not None:
    coverage_co = args.cov
else:
    coverage_co = 0.0
if args.cval is not None:
    cval_co = args.cval
else:
    cval_co = 0.01
if args.ival is not None:
    ival_co = args.ival
else:
    ival_co = 0.01
if args.acc is not None:
    acc_co = args.acc
else:
    acc_co = 0.6 


#==============================================================================
#  IO
#==============================================================================
proteome_filename = args.seqdb #'/home/nicolas.chevrollier/herbiFun_project/hmmbuilder/datas/mgg_70-15_8.fasta'
input_msa = args.fasta #'/home/nicolas.chevrollier/herbiFun_project/hmmbuilder/datas/KS.msa'


sepdir = os.path.dirname(input_msa)
print sepdir
if len(sepdir) != 0:
    dom_name = input_msa.split(sepdir)[-1].split('/')[-1].split('.')[0]
else:
    dom_name = input_msa.split('.')[0]
outname = dom_name+'_hmmbuild'


path2output = os.path.abspath(args.dir)+'/'+outname+'_'+date+'/'


if not os.path.exists(path2output):
    os.makedirs(path2output)

runs_output = 'runs_output/'


#==============================================================================
# Log file handler
#==============================================================================
logger = logger.define_logger(path2output=path2output, outname=outname, date=date)


#==============================================================================
# Inputs report
#==============================================================================
working_dir = os.getcwd()
cmd_input = 'python '+' '.join(sys.argv)+'\n'

txt = '{}'.format(40*'-')
txt = txt + ' Starting hmmbuilder '.format()
txt = txt + '{}\n\n'.format(40*'-')
logger.info(txt+cmd_input)


#==============================================================================
# Main
#==============================================================================
is_newdomain = 1
seqnb_tmp = seqio.seqnumber(input_msa)
iter_nb = 1
n_conv = 0
while is_newdomain:
    if iter_nb == 1:
        clust=True
    else:
        clust=False
    
    # Building hmm        
    msa, clw, hmm, hmmsearch_dic = buildhmm(logger=logger, input_msa=input_msa, proteome_filename=proteome_filename, 
                                                    path2output=path2output+runs_output, iter_nb=iter_nb, clust=clust,
                                                    identity=identity, coverage_co=coverage_co, cval_co=cval_co, ival_co=ival_co, acc_co=acc_co)
    
    seqnb = seqio.seqnumber(msa)
    
    txt = '\n{} hits found and accepted ({} in the previous iteration)\n'.format(seqnb, seqnb_tmp) 
    logger.info(txt)
    
    if seqnb != seqnb_tmp and n_conv <= 3:
        seqnb_tmp = seqnb
        iter_nb += 1
        
        if abs(seqnb-seqnb_tmp) == 1:
            n_conv+=1           
    else:
        is_newdomain = 0
        if n_conv > 3:
            txt = '{}'.format(40*'-')
            txt = txt + ' Computation has finished ! The convergence criteria has been reached. '.format()
            txt = txt + '{}'.format(40*'-')
            logger.info(txt)
        else:
            txt = '{}'.format(40*'-')
            txt = txt + ' Computation has finished ! No new hits have been found '.format()
            txt = txt + '{}'.format(40*'-')
            logger.info(txt)


#==============================================================================
# final ouptuts
#==============================================================================
msa_prev = msa.replace(path2output+runs_output+dom_name+'-'+str(iter_nb+1), path2output+runs_output+dom_name+'-'+str(iter_nb))

hmm_final = path2output+dom_name+'.hmm'
clw_final = path2output+dom_name+'.seed'
msa_final = path2output+dom_name+'.msa'

shutil.move(hmm, hmm_final)
shutil.move(clw, clw_final)
shutil.move(msa_prev, msa_final)

txt = '\nOutput files:\n'.format()
txt = txt + '\t- HMM profile: {}\n\t- Seeds\'alignment: {}\n\t- Seed sequences: {}\n'.format(hmm_final, clw_final, msa_final)
logger.info(txt)






























