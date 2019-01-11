# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 18:19:10 2018

@author: nicolas.chevrollier
"""

def is_domrelevant(domtbl, coverage_co=0.0, cval_co=0.01, ival_co=0.01, acc_co=0.6):
    """
    Input: domtbl is an instance of the Domtblparser object
    Output: boolean (1 if domain meets criteria, 0 otherwise)
    """
    is_relevant = 0 
    
    envlen = domtbl.env_to - domtbl.env_from + 1
    if float(envlen) >= domtbl.qlen * coverage_co:
        if domtbl.dom_cval <= cval_co and domtbl.dom_ival <= ival_co and domtbl.acc >= acc_co:
            is_relevant = 1
            
    return is_relevant


def filtering_hits(domtbl=None, coverage_co=0.0, cval_co=0.01, ival_co=0.01, acc_co=0.6, logger=None):
    """
    Remove domains that don't pass criteria selection, and hits if no domain persists
    """
    for hit in sorted(domtbl):
        for dom_id in sorted(domtbl[hit]):
            if is_domrelevant(domtbl[hit][dom_id], coverage_co=coverage_co, cval_co=cval_co, ival_co=ival_co, acc_co=acc_co):        
                if logger is None: pass
#                    print('ACCEPTED -> {}'.format(domtbl[hit][dom_id].printdom()))
                else: 
                    logger.info('ACCEPTED -> {}'.format(domtbl[hit][dom_id].printdom()))
            else:
                if logger is None: pass
#                    print('REJECTED -> {}'.format(domtbl[hit][dom_id].printdom()))
                else:
                    logger.info('REJECTED -> {}'.format(domtbl[hit][dom_id].printdom()))
                del domtbl[hit][dom_id]
        if len(domtbl[hit]) == 0:
            del domtbl[hit]


class Domtblparser(object):
    def __init__(self, line, domtype='search'):        
        cols = line.split()[:22]
        if domtype == 'search':
            self.tname = cols[0]
            self.tlen = int(cols[2])
            self.qname = cols[3]
            self.qlen = int(cols[5])
        elif domtype == 'scan':
            self.tname = cols[3]
            self.tlen = int(cols[5])
            self.qname = cols[0]
            self.qlen = int(cols[2])
        self.seq_eval = float(cols[6])
        self.seq_score = float(cols[7])
        self.seq_bias = float(cols[8])
        self.dom_id = int(cols[9])
        self.dom_nb = int(cols[10])
        self.dom_cval = float(cols[11])
        self.dom_ival = float(cols[12])
        self.dom_score = float(cols[13])
        self.dom_bias = float(cols[14])
        self.hmm_from = int(cols[15])
        self.hmm_to = int(cols[16])
        self.ali_from = int(cols[17])
        self.ali_to = int(cols[18])
        self.env_from = int(cols[19])
        self.env_to = int(cols[20])
        self.acc = float(cols[21])
        self.sequence = None
        
    def printdom(self):
        return 'tname: {}\tqname: {}\tid: {}\tcval: {}\tival: {}\tscore: {}\tacc: {}\tenvlen: {}\talilen: {}\thmmlen: {}\tqlen: {}'.format(self.tname, self.qname, self.dom_id,
              self.dom_cval, self.dom_ival, self.dom_score, self.acc, self.get_envlen(), self.get_alilen(), self.get_hmmlen(), self.qlen)
        
    def get_envlen(self):
        return self.env_to-self.env_from+1
        
    def get_alilen(self):
        return self.ali_to-self.ali_from+1
        
    def get_hmmlen(self):
        return self.hmm_to-self.hmm_from+1

        
def domtbl_dic(domtblout=None, domtype='search'):
    """
    Input: domtblout hmmer format file
    Output: dictionary of dictionary
               - dic[target_name][domain_id] = Domtblparser instance
    """
    if domtblout is not None:
        hmm_hits = {}
        with open(domtblout, "r") as hmmsearch_file:
            for l in hmmsearch_file:
                if not l.startswith("#"):
                    domtbl = Domtblparser(l, domtype)
                    tname = domtbl.tname
                    key = domtbl.qname+'_'+str(domtbl.dom_id)

                    if tname not in hmm_hits:
                        hmm_hits[tname] = {}
                    hmm_hits[tname][key] = domtbl                  
                    
        return hmm_hits
    else: 
        print('A domtblout hmmer format is required...'.format())

  
def get_domtblseq(domtbl_dic=None, fasta=None):
    """
    Inputs: - output of domtbl_dic()
            - output of seqio.read_fasta() (dictionary)
    Function: Gets sequences of hmmsearch hits
    """
    for hit_id in sorted(domtbl_dic):
        for dom_id in domtbl_dic[hit_id]:
            coor_from = domtbl_dic[hit_id][dom_id].ali_from - 1
            coor_to = domtbl_dic[hit_id][dom_id].ali_to
            dom_seq = fasta[hit_id][coor_from:coor_to]            
            
            domtbl_dic[hit_id][dom_id].sequence = dom_seq


def get_domtblfasta(hmmsearch_hits):
    """
    Input: output of domtbl_dic()
    Output: list of fasta sequences' hits
    """
    dom_fasta = []
    for hit_id in sorted(hmmsearch_hits):
        for i,dom_id in enumerate(hmmsearch_hits[hit_id]):
            header = '>{}_{}\n'.format(hit_id, i+1)
            sequence =  '{}\n'.format(hmmsearch_hits[hit_id][dom_id].sequence)
            dom_fasta.append(header+sequence)
    
    return dom_fasta



if __name__ == '__main__':
    
    scanoutput = '/home/nicolas.chevrollier/herbiFun_project/hmmbuilder/datas/mgg_70-15_8.domtblout'
    
    # Gets hmm scan hits
    hmmscan_hits = domtbl_dic(domtblout=scanoutput, domtype='scan')
    # Keeps only domains passing the criteria
    filtering_hits(domtbl=hmmscan_hits, coverage_co=0.0, cval_co=0.01, ival_co=0.01, acc_co=0.6, logger=None)























