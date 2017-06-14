#!/usr/bin/env python3

import sys
import re
from operator import itemgetter

class TrnAllele:
    def __init__(self, allele, prefix, suffix, mtx, trn):
        self.id = allele.id
        self.allele = allele
        self.prefix = prefix
        self.suffix = suffix
        self.context = prefix.id + ':' + (suffix.id if suffix.id != '-' else '')
        self.mtx = mtx
        self.trn = trn
# end of TrnAllele

class TrnVariation:
    def __init__(self, snp, beg, end):
        self.id = snp.id
        self.ref_allele_id = None
        self.snp = snp
        self.beg_vertebra = beg
        self.end_vertebra = end
        self._alleles = {}
        self._alleles_in_context = {}
        self._nonsyn_alleles = {}
        self._allele_usage = {}
        self._allele_rank = {}
    
    def appendTrnAllele(self, trn_allele):
        if not trn_allele.context in self._alleles_in_context:
            self._alleles_in_context[trn_allele.context] = {}
        if not trn_allele.trn in self._alleles_in_context[trn_allele.context]:
            self._alleles_in_context[trn_allele.context][trn_allele.trn] = []
        self._alleles_in_context[trn_allele.context][trn_allele.trn].append(trn_allele)
        if trn_allele.allele.isReference():
            self.ref_allele_id = trn_allele.id
        if not trn_allele.id in self._alleles:
            self._alleles[trn_allele.id] = []
        self._alleles[trn_allele.id].append(trn_allele)
        
    def getTrnAlleles(self):
        return self._alleles
    
    def getNonSynAlleles(self):
        return self._nonsyn_alleles
    
    def joinSynonymAlleles(self):
        context_list = list(self._alleles_in_context.keys())
        context_list.sort(key=len)
        context2remove = set()
        for i in range(len(context_list)):
            if i + 1 < len(context_list):
                for j in range(i + 1, len(context_list)):
                    if re.match(context_list[i], context_list[j]):
                        for trn in self._alleles_in_context[context_list[i]]:
                            if trn not in self._alleles_in_context[context_list[j]]:
                                self._alleles_in_context[context_list[j]][trn] = []
                            self._alleles_in_context[context_list[j]][trn].extend(self._alleles_in_context[context_list[i]][trn])
                        context2remove.add(context_list[i])
        for context in context2remove:
            if len(self._alleles_in_context[context]) < 2:
                self._alleles_in_context.pop(context)
        
        n = 0
        for context in self._alleles_in_context:
            for trn in self._alleles_in_context[context]:
                n += 1
                l = len(self._alleles_in_context[context][trn])
                if l < 1:
                    sys.exit("ERROR: wrong allele context - {}\n".format(self.id))
                l = 1 / l
                for trn_allele in self._alleles_in_context[context][trn]:
                    if trn_allele.id in self._allele_usage:
                        self._allele_usage[trn_allele.id] += l
                    else:
                        self._allele_usage[trn_allele.id] = l
        if n:
            #sys.exit("ERROR: wrong translated variation - {}\n".format(self.id))
            allele_ranks = []
            for allele_id in self._allele_usage:
                self._allele_usage[allele_id] /= n # allele usage calculation
                allele_ranks.append({'id': allele_id, 'usage': self._allele_usage[allele_id], 'ref': 1 if allele_id == self.ref_allele_id else 0})
            allele_ranks = sorted(allele_ranks, key=itemgetter('usage', 'ref'))
            for i in range(len(allele_ranks)):
                self._allele_rank[allele_ranks[i]['id']] = i
            
            for context in self._alleles_in_context:
                for trn in self._alleles_in_context[context]:
                    max_rank = -1
                    best_allele_id = None
                    prefix_seq = None
                    for trn_allele in self._alleles_in_context[context][trn]:
                        if max_rank < self._allele_rank[trn_allele.id]:
                            max_rank = self._allele_rank[trn_allele.id]
                            best_allele_id = trn_allele.id
                            prefix_seq = trn_allele.prefix.seq
                    if not best_allele_id:
                        sys.exit("ERROR: no best allele - {}\n".format(self.id))
                    if len(self._alleles_in_context[context]) > 1:
                        self._nonsyn_alleles[best_allele_id] = self._alleles[best_allele_id]
                        self._nonsyn_alleles[best_allele_id][0].allele.setNonSyn()
                    for node in self.beg_vertebra.getNodes():
                        if node.pos == 152312604:
                            ddd=1
                        for allele in node.getAlleles():
                            if allele.isPhased() or allele.id == best_allele_id:
                                node.appendPrefixSeq(prefix_seq)
                                break
        return self._nonsyn_alleles
    
    def isNonSynonymous(self, allele_id = None):
        if allele_id:
            if allele_id in self._nonsyn_alleles:
                return True
            return False
        if len(self._nonsyn_alleles) > 1:
            return True
        return False
# end of TrnVariation
