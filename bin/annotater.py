#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 17:33:55 2018

@author: nicolas.chevrollier
"""
#==============================================================================
# Imports
#==============================================================================
import os

#import _path2lib
#_path2lib.link_bin2lib()

#import sys
#sys.path.append('/home/nicolas.chevrollier/herbiFun_project/hmmbuilder')

import lib.external_prog as external_prog
import lib.seqio as seqio
import lib.hmmdomtbl as hmmdomtbl

from lxml import etree
import argparse


#==============================================================================
# Functions
#==============================================================================
def annotation_rules(rule_filename=None):
    rule_dic = {'Rules': {}}
    with open(rule_filename, 'r') as rule_file:
        for l in rule_file:
            if l.startswith('#'):
                continue
            classtype = l.split('|')[0].strip()
            name = l.split('|')[1].strip()
            rule = l.split('|')[2].strip()
            
            rule_dic['Rules'][rule] = {'Class': classtype, 'Name': name}
    
    return rule_dic

def get_xml(annotated_proteins=None, hmmscan_hits=None, proteome_name=None, proteome_fasta=None):
    """
    - Returns a dictionary: xml_prot[hit] = xml
    """
    xml_prot = {}
    for hit in sorted(annotated_proteins):   
        protein = etree.Element("protein")
        protein_id = etree.SubElement(protein, "proteome")
        protein_id.text = os.path.basename(proteome_name)
        protein_id = etree.SubElement(protein, "id")
        protein_id.text = hit
        protein_class = etree.SubElement(protein, "class")
        protein_class.text = annotated_proteins[hit]
        protein_seq = etree.SubElement(protein, "sequence")
        protein_seq.text = proteome_fasta[hit]
        protein_length =  etree.SubElement(protein, "length")
        protein_length.text = str(len(proteome_fasta[hit]))
        
        for i in sorted([ (hmmscan_hits[hit][x].qname, hmmscan_hits[hit][x].ali_from, 
                    hmmscan_hits[hit][x].ali_to, hmmscan_hits[hit][x].dom_score, hmmscan_hits[hit][x].dom_ival) for x in hmmscan_hits[hit] ],  key=lambda x: x[1]):
            name =  i[0]
            coor_from = i[1]
            coor_to = i[2]
            length = coor_to - coor_from + 1
            score = i[3]
            ival = i[4]
            
            domain = etree.SubElement(protein, "domain")    
            domain_name = etree.SubElement(domain, "name")
            domain_name.text = name
            domain_length = etree.SubElement(domain, "length")
            domain_length.text = str(length)
            domain_from = etree.SubElement(domain, "from")
            domain_from.text = str(coor_from)
            domain_to = etree.SubElement(domain, "to")
            domain_to.text = str(coor_to)
            domain_score = etree.SubElement(domain, "score")
            domain_score.text = str(score)
            domain_ival = etree.SubElement(domain, "i-eval")
            domain_ival.text = str(ival)    
    
        xml = (etree.tostring(protein, pretty_print=True))
        xml_prot[hit] = xml
        
    return xml_prot


def get_annotatedprot(hmmscan_hits=None, rule_dic=None):
    """
    - Annotates proteins according to the rules
    - Returns a dictionary: annotated_prot[hit] = prot_class
    """
    
    annotated_prot = {}
    for hit in sorted(hmmscan_hits):
        prot_class = 'NA'
        dom_list = [ hmmscan_hits[hit][x].qname for x in hmmscan_hits[hit] ]
        for rule in sorted(rule_dic['Rules'], key=lambda x: len(x), reverse=True):
            rule_list = [ x.strip() for x in rule.split(',') ]
            if len(set(dom_list) & set(rule_list)) == len(rule_list):            
                prot_class = rule_dic['Rules'][rule]['Class']
                break
        if prot_class is not 'NA':
            annotated_prot[hit] = prot_class
            #print hit, prot_class, [ (hmmscan_hits[hit][x].qname, hmmscan_hits[hit][x].dom_score) for x in hmmscan_hits[hit] ]
            
    return annotated_prot


def write_xml(path2output=None, xml_dic=None, annotated_proteins=None):    
    for hit in sorted(xml_dic):
        prot_class = annotated_proteins[hit]
        xml_filename = path2output+'/'+prot_class+'_'+hit+'.xml'
        
        with open(xml_filename, 'w') as xml_file:
            xml_file.write(xml_dic[hit])


#==============================================================================
#  Argument parsing
#==============================================================================
parser = argparse.ArgumentParser(description='Iterative building of hmm profiles')
parser.add_argument("-proteome", required = True, nargs = "?", help = "Proteome fasta file")
parser.add_argument("-hmmdb", required = True, nargs = "?", help = "HMM profile database")
parser.add_argument("-rules", required = True, nargs = "?", help = "File containing rules")
parser.add_argument("-dir", required = True, nargs = "?", help = "Output directory")
parser.add_argument("-cov", required = False, help = "Minimum percentage of coverage alignment between hmm hit and hmm profile (0.0)", type=float)
parser.add_argument("-cval", required = False, help = "hmmer conditional e-value cutoff (0.01)", type=float)
parser.add_argument("-ival", required = False, help = "hmmer independant e-value cutoff (0.01)", type=float)
parser.add_argument("-acc", required = False, help = "hmmer mean probability of the alignment accuracy between each residues of the target and the corresponding hmm state (0.6)", type=float)
args = parser.parse_args()

# Optional parameters
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

proteome_name = args.proteome
hmmdatabase = args.hmmdb
xml_outputs = args.dir
rule_filename = args.rules

#proteome_name = '/home/nicolas.chevrollier/herbiFun_project/hmmbuilder_v2/datas/mgg_70-15_8.fasta'
#hmmdatabase = '/home/nicolas.chevrollier/herbiFun_project/hmmbuilder_v2/datas/database.hmm'
#scanoutput = '/home/nicolas.chevrollier/herbiFun_project/hmmbuilder_v2/datas/mgg_70-15_8.domtblout'
#xml_outputs = '/home/nicolas.chevrollier/herbiFun_project/hmmbuilder_v2/datas/annotated_proteins/'
#rule_filename = '/home/nicolas.chevrollier/herbiFun_project/hmmbuilder_v2/datas/annotation.rules'


scanoutput = proteome_name.split('.fasta')[0]+'.domtblout'
external_prog.hmmscan(proteome_filename=proteome_name, scanoutput=scanoutput, hmmdb_filename=hmmdatabase)

# Gets hmm scan hits
hmmscan_hits = hmmdomtbl.domtbl_dic(domtblout=scanoutput, domtype='scan')
# Keeps only domains passing the criteria
hmmdomtbl.filtering_hits(domtbl=hmmscan_hits, coverage_co=coverage_co, cval_co=cval_co, ival_co=ival_co, acc_co=acc_co, logger=None)

       
# Gets annotation rules
rule_dic = annotation_rules(rule_filename=rule_filename)

#annotated_prot = {}
#for hit in sorted(hmmscan_hits):
#    prot_class = 'NA'
#    dom_list = [ hmmscan_hits[hit][x].qname for x in hmmscan_hits[hit] ]
#    for rule in sorted(rule_dic['Rules'], key=lambda x: len(x), reverse=True):
#        rule_list = [ x.strip() for x in rule.split(',') ]
#        if len(set(dom_list)) == len(rule_list):
#            if len(set(dom_list) & set(rule_list)) == len(rule_list):            
#                prot_class = rule_dic['Rules'][rule]['Class']
#                break
#        else: pass
#    if prot_class is not 'NA':
#        annotated_prot[hit] = prot_class
#
#rule_list = [ x.strip() for x in 'PP'.split(',') ]
#dom_list = ['PP', 'KS']
#
#for hit in annotated_prot:
#    if annotated_prot[hit] == 'dom_PP':
#        print hit, [ hmmscan_hits[hit][x].qname for x in hmmscan_hits[hit] ]


# Annotates proteins according to the rules
annotated_proteins = get_annotatedprot(hmmscan_hits=hmmscan_hits, rule_dic=rule_dic)        

# Gets sequences of domains
proteome_fasta = seqio.read_fasta(sequences=proteome_name)
hmmdomtbl.get_domtblseq(domtbl_dic=hmmscan_hits, fasta=proteome_fasta) 

# Formats annotated proteins into xml format
xml_prot = get_xml(annotated_proteins=annotated_proteins, hmmscan_hits=hmmscan_hits, proteome_name=proteome_name, proteome_fasta=proteome_fasta)

# Writes xml files
if not os.path.exists(xml_outputs):
    os.makedirs(xml_outputs)

write_xml(path2output=xml_outputs, xml_dic=xml_prot, annotated_proteins=annotated_proteins)




























