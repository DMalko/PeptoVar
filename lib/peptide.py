#!/usr/bin/env python3

import math
import sys
from copy import copy
from lib.seqtools import translate, is_stopcodon

DEBUG = 1

class Peptide:
    def __init__(self, chrom, transcript_id, sample, fshifts, mtx, trn, beg, end, alleles):
        self.chrom = chrom
        self.trn_id = transcript_id
        self.sample = sample
        self.fshifts = fshifts
        self.mtx = mtx
        self.trn = trn
        self.beg = beg
        self.end = end
        self._alleles = alleles
    
    def getTranslation(self, aa_length = None):
        if not len(self.matrix) or len(self.matrix) % 3:
            sys.exit("ERROR: wrong peptide length - {}\n".format(self.matrix))
        return (translate(self.matrix))[0]
    
    def appendAllele(self, allele, pos):
        if allele.id not in self._alleles:
            self._alleles[allele.id] = {}
            self._alleles[allele.id]['beg'] = pos
            self._alleles[allele.id]['end'] = pos
        else:
            self._alleles[allele.id]['end'] = pos
    
    def getAlleles(self):
        return self._alleles
    
    def _strpos(self, pos):
        if pos:
            if pos['overpos']:
                return "{}+{}({})".format(pos['pos'], pos['overpos'], pos['id'])
            else:
                return str(pos['pos'])
        return '0'
    
    def getBegPos(self):
        return self._strpos(self.beg)
    
    def getEndPos(self):
        return self._strpos(self.end)
    
    def format(self):
        fsh = []
        for fshift in self.fshifts:
            fsh.append(fshift.id)
        var = []
        for allele_id, pos in sorted(self.getAlleles().items(), key=lambda x: x[1]['beg']):
            allele_str = '({}..{}){}'.format(pos['beg'], pos['end'], allele_id)
            if pos['mark_as_synonymous']:
                allele_str = '[' + allele_str + ']'
            var.append(allele_str)
        if not len(var):
            var.append('-')
        
        return [self.chrom, self.trn_id, self.sample.name, self.sample.allele_1, self.sample.allele_2, self.getBegPos(), self.getEndPos(), ';'.join(fsh), ';'.join(var), self.trn, self.mtx]
# end of Peptide

class PeptideDriver:
    def __init__(self, chrom, transcript_id, sample, length_set, opt = False):
        self._chrom = chrom
        self._trn_id = transcript_id
        self._sample = sample
        self._start_node = None
        self._end_node = None
        self._len_set = sorted(length_set, reverse=True)
        self._len = self._len_set.pop() if len(self._len_set) else None
        self._mtx = ''
        self._beg = None
        self._end = None
        self._loc_pos = 0
        self._alleles = {}
        self._peptides = {}
        self._optimization = opt
    
    def copy(self):
        new_driver = copy(self)
        self._alleles = self._alleles.copy()
        return new_driver
    
    def setPeptContainer(self, container):
        self._peptides = container
    
    def getPeptContainer(self):
        return self._peptides
    
    def getNext(self):
        if self._end_node:
            prefix = None
            if self._optimization:
                prefix_len = len(self._mtx) % 3
                prefix = self._mtx[-prefix_len:] if prefix_len else ''
            return self._end_node.getNext(self._sample, prefix)
        else:
            return []
    
    def appendNode(self, node):
        if self._len is not None:
            pos = node.getPos()
            # peptide position is determined as a "variation shadow" on the genome sequence
            if self._loc_pos == 0:
                self._start_node = node
                self._beg = pos
                self._end = pos.copy()
            else:
                if pos['overpos']:
                    if self._end['overpos']:
                        self._end = pos
                    else:
                        self._end['pos'] += 1
                else:
                    if self._beg['overpos']:
                        self._beg['pos'] = node.pos - self._loc_pos
                        self._beg['overpos'] = pos['overpos']
                        self._beg['id'] = pos['id']
                    self._end = pos
            self._end_node = node
            self._loc_pos += 1
            
            self._mtx += node.nucl
            mtx_len = len(self._mtx)
            trn_len = mtx_len / 3
            
            for allele in node.getAlleles(self._sample):
                if allele.id not in self._alleles:
                    self._alleles[allele.id] = {}
                    self._alleles[allele.id]['beg'] = mtx_len
                    self._alleles[allele.id]['end'] = mtx_len
                    if self._optimization and not allele.isNonSyn():
                        self._alleles[allele.id]['mark_as_synonymous'] = True
                    else:
                        self._alleles[allele.id]['mark_as_synonymous'] = False
                else:
                    self._alleles[allele.id]['end'] = mtx_len
            
            get_stop = 0
            if mtx_len > 2 and mtx_len % 3 == 0 and is_stopcodon(self._mtx[-3:]):
                get_stop = 1
            
            if self._len != 0 and get_stop:
                return False
            
            if self._len != 0 and trn_len == self._len or self._len == 0 and (get_stop or len(self.getNext()) == 0):
                peptide = Peptide(self._chrom, self._trn_id, self._sample, self._start_node.getFShiftPathSet(self._sample), self._mtx, (translate(self._mtx))[0], self._beg, self._end, self._alleles.copy())
                pept_record = peptide.format()
                if DEBUG:
                    print(" ".join(str(item) for item in pept_record))
                
                if self._sample.name not in self._peptides:
                    self._peptides[self._sample.name] = []
                self._peptides[self._sample.name].append(pept_record)
                
                if len(self._len_set):
                    self._len = self._len_set.pop()
                    return True
                else:
                    self._len = None
                    return False
            elif self._len > 0 and trn_len < self._len or self._len == 0:
                return True
            else:
                sys.exit("ERROR: wrong peptide length - {}\n".format(self._mtx))
        return False
# end of PeptideContainer
