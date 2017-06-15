#!/usr/bin/env python3

import re
from lib.seqtools import translate, complement

class Seq:
    def __init__(self, strand, beg, end, seq):
        self._seq = list(seq)
        self.strand = strand
        self.beg = beg
        self.end = end
    
    def setSeq(self, lbeg, lend, seq):
        llen = lend - lbeg + 1
        if llen > 0:
            seq = list(seq)
            if len(seq) < llen:
                seq.extend('-' for i in range (llen - len(seq)))
            elif len(seq) > llen:
                tail = seq[llen - 1:]
                seq = seq[:llen - 1]
                seq.append(tail)
            if len(seq) != llen:
                print("ERROR: wrong allele_seq length {}..{}".format(lbeg, lend))
                exit()
        else:
            print("ERROR: seq.beg > seq.end {}..{}".format(lbeg, lend))
            exit()
        
        for i in range(len(seq)):
            if self._seq[lbeg + i] != '-': # some VCF files have wrong overlapped alleles
                self._seq[lbeg + i] = seq[i]
    
    def getSeq(self, beg = 0, end = None):
        seq = ''
        if end is None:
            end = len(self._seq)
        else:
            end += 1
        for nucl in self._seq[beg:end]:
            if type(nucl) is list:
                seq += ''.join(nucl)
            elif nucl != '-':
                seq += nucl
        return seq
# end of Seq

class Sequence:
    def __init__(self, chrom, trn_id, samples):
        self.chrom = chrom
        self.trn_id = trn_id
        self.strand = None
        self._smplseq = {}
        self._refseq = []
        self._smplprot = {}
        self._smplpept = {}
        self.refprot = ''
        self.samples = samples
        self._strangers = {}
    
    def append(self, strand, beg, end, seq):
        if not self.strand:
            self.strand = strand
        elif self.strand != strand:
            print("ERROR: Wrong strand!")
        for sample in self.samples:
            if sample.id not in self._smplseq:
                self._smplseq[sample.id] = []
            self._smplseq[sample.id].append(Seq(strand, beg, end, seq))
            self._refseq.append(Seq(strand, beg, end, seq))
    
    def modify(self, snp):
        for sample in self.samples:
            for allele in snp.getSampleAlleles(sample):
                if allele.isReference():
                    for exseq in self._refseq:
                        if exseq.beg <= snp.beg and snp.end <= exseq.end:
                            lbeg = snp.beg - exseq.beg
                            lend = snp.end - exseq.beg
                            refseq = exseq.getSeq(lbeg, lend)
                            if refseq.upper() != allele.seq:
                                print("Oy Gevalt! Wrong reference allale!")
                            break
                else:
                    for exseq in self._smplseq[sample.id]:
                        if exseq.beg <= snp.beg and snp.end <= exseq.end:
                            lbeg = snp.beg - exseq.beg
                            lend = snp.end - exseq.beg
                            exseq.setSeq(lbeg, lend, allele.seq)
                            break
    
    def translate(self):
        self._refseq = sorted(self._refseq, key = lambda x: x.beg)
        seq = "".join(x.getSeq() for x in self._refseq)
        if self.strand == '-':
            self.refprot = (translate(complement(seq[::-1])))[0]
        else:
            self.refprot = (translate(seq))[0]
        
        for sample in self.samples:
            self._smplseq[sample.id] = sorted(self._smplseq[sample.id], key = lambda x: x.beg)
            
            seq = "".join(x.getSeq() for x in self._smplseq[sample.id])
            if self.strand == '-':
                protein = (translate(complement(seq[::-1])))[0]
                protein = re.sub(r"\*.*", "", protein)
                self._smplprot[sample.id] = protein
            else:
                protein = (translate(seq))[0]
                protein = re.sub(r"\*.*", "", protein)
                self._smplprot[sample.id] = protein
        return self._smplprot
    
    def peptides(self, length_set):
        for sample_id in self._smplprot:
            for length in length_set:
                if length: # length == 0 - full length protein
                    for i in range(len(self._smplprot[sample_id])):
                        pept = self._smplprot[sample_id][i : i + length]
                        if len(pept) == length:
                            if sample_id not in self._smplpept:
                                self._smplpept[sample_id] = set()
                            self._smplpept[sample_id].add(pept)
                        else:
                            break
    
    def refpeptides(self, length):
        peptset = set()
        for i in range(len(self.refprot)):
            pept = self.refprot[i : i + length - 1]
            if len(pept) == length:
                peptset.add(pept)
            else:
                break
        return peptset
    
    def add2compare(self, sample_peptides):
        for item in sample_peptides:
            sample_id = "_".join((item[2], str(item[3]), str(item[4])))
            peptide = item[9]
            if sample_id not in self._strangers:
                self._strangers[sample_id] = []
            self._strangers[sample_id].append(peptide)
    
    def compare(self):
        misspept = []
        for sample_id, peptides in self._strangers.items():
            for peptide in peptides:
                if sample_id not in self._smplpept or peptide not in self._smplpept[sample_id]:
                    misspept.append("graphpept\t{}\t{}\t{}\t{}\n".format(self.chrom, self.trn_id, sample_id, peptide))
        for sample_id, peptides in self._smplpept.items():
            for peptide in peptides:
                if sample_id not in self._strangers or peptide not in self._strangers[sample_id]:
                    misspept.append("dummypept\t{}\t{}\t{}\t{}\n".format(self.chrom, self.trn_id, sample_id, peptide))
        return misspept
# end of Sequence
