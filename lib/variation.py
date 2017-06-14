#!/usr/bin/env python3

import sys
from lib.seqtools import complement
from lib.sample import SampleContainer, Sample

class Allele:
    def __init__(self, snp_id, seq, strand, freq, ref = None):
        self.snp_id = snp_id
        self.seq = seq.upper()
        self.length = len(seq)
        self.freq = freq
        self.strand = strand
        self.id = ":".join((self.snp_id.lower(), self.seq.upper())) + ("(ref)" if ref else "(alt)")
        self._is_ref = True if ref else False
        self._is_nonsyn = False
        self._indel = False
        self._samples = SampleContainer()
        self._phased = SampleContainer() # samples for which the allele is phased
    
    def addSample(self, sample, phased = False):
        if phased:
            self._phased.appendSample(sample)
        return self._samples.appendSample(sample)
    
    def getSamples(self):
        return self._samples.getSample()
    
    def getSample(self, sample):
        return self._samples.getSample(sample)
    
    def isPhased(self, sample = None):
        if sample:
            return self._phased.getSample(sample)
        return len(self._phased.getSample())
    
    def isReference(self):
        return self._is_ref
    
    def setInDel(self):
        self._indel = True
        
    def setNonSyn(self):
        self._is_nonsyn = True
    
    def isNonSyn(self):
        return self._is_nonsyn
    
    def isInDel(self):
        return self._indel
        
# end of Allele

class SNP:
    def __init__(self, snp, samples):
        if hasattr(snp, 'info'): # normalization
            if 'MATCHED_REV' in snp.info:
                if snp.info['MATCHED_REV'] == True:
                    new_alt = []
                    for allele in snp.alleles:
                        new_alt.append(complement(allele)[::-1])
                    snp.alleles = tuple(new_alt)
        
        self.id = snp.id
        self.beg = snp.start + 1
        self.end = snp.stop
        self.strand = '+'
        self.ref = Allele(snp.id, snp.ref, self.strand, 0, 'reference')
        self.alleles = [self.ref] # all alleles in the variation (self.alleles[0] - the reference allele forever!!!)
        self._sample_allele_index = set() # self.allele indexes for samples (they will be taken for graph construction)
        self._indel = False # insertion/deletion variation type
        
        ref_af = -1
        if hasattr(snp, 'info'):
            if 'AF' in snp.info: # according to VCF file description allele frequency (AF) specified ONLY for alternative alleles
                if len(snp.alts) == len(snp.info['AF']):
                    ref_af = 1
                    alt_alleles = []
                    for i, allele in enumerate(snp.alts):
                        if allele == snp.ref: continue
                        af = round(float(snp.info['AF'][i]), 4)
                        ref_af -= af # AF for reference allele should be calculated
                        new_allele = Allele(snp.id, allele, self.strand, af)
                        indel = (len(self.ref.seq) - len(allele)) % 3
                        if indel:
                            self._indel = True
                            new_allele.setInDel()
                        alt_alleles.append(new_allele)
                    if ref_af < 0:
                        print('WARNING: wrong reference allele frequency - SNP {}'.format(snp.id))
                    else:
                        self.ref.freq = ref_af
                        self.alleles.extend(alt_alleles)
                else:
                    print('WARNING: wrong allele frequency data in AF field - SNP {}'.format(snp.id))
        if ref_af < 0:
            # set AF = 0 for all alleles
            for allele in snp.alts:
                # self.ref has been initialized by '0'
                new_allele = Allele(snp.id, allele, self.strand, 0)
                frame_shift = (len(self.ref.seq) - len(allele)) % 3
                if frame_shift:
                    self._indel = True
                    new_allele.setFrameShift()
                self.alleles.append(new_allele)
        
        for sample in samples:
            if sample == 'virtual': # the virtual sample (all alleles are unphased)
                for i, allele in enumerate(self.alleles):
                    allele.addSample(Sample(sample, 1, 1))
                    self._sample_allele_index.add(i)
            else:
                if hasattr(snp, 'samples'):
                    try:
                        sample in snp.samples
                    except:
                        raise
                    else:
                        allele1_idx = self._findAlleleIndex(snp.samples[sample].alleles[0])
                        allele2_idx = self._findAlleleIndex(snp.samples[sample].alleles[1])
                        if allele1_idx is None or allele2_idx is None:
                            sys.exit("ERROR: can not find the snp allele ({}) for the sample {}\n".format(self.id, sample))
                        if hasattr(snp.samples[sample], 'phased') and snp.samples[sample].phased:
                            self.alleles[allele1_idx].addSample(Sample(sample, 1, 0), True) # True is phased
                            self.alleles[allele2_idx].addSample(Sample(sample, 0, 1), True) # True is phased
                        else:
                            self.alleles[allele1_idx].addSample(Sample(sample, 1, 1))
                            self.alleles[allele2_idx].addSample(Sample(sample, 1, 1))
                        self._sample_allele_index.update([allele1_idx, allele2_idx])
                else:
                    raise ValueError("No any samples in the VCF file\n")
    
    def _findAlleleIndex(self, seq): # find the allele with a given sequence
        for i, allele in enumerate(self.alleles):
            if allele.seq == seq:
                return i
        return None
    
    def setComplement(self):
        self.strand = '-' if self.strand == '+' else '+'
        for allele in self.alleles:
            allele.seq = complement(allele.seq)
            allele.strand = self.strand
    
    def filterAF(self, min_af):
        if min_af > 0:
            self._sample_allele_index = set(i for i in self._sample_allele_index if self.alleles[i].freq >= min_af)
            if len(self._sample_allele_index) > 1:
                return True
            elif len(self._sample_allele_index) == 1 and self._sample_allele_index[0] != 0: # reference allele has index = 0
                return True
            else:
                return False
        return True
    
    def getSampleAlleles(self, sample = None):
        alleles = []
        for i in self._sample_allele_index:
            if sample:
                if self.alleles[i].getSample(sample):
                    alleles.append(self.alleles[i])
            else:
                alleles.append(self.alleles[i])
        return alleles
    
    def isInDel(self):
        return self._indel
# end of SNP
