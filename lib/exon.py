#!/usr/bin/env python3

import sys
import weakref
from lib.graph import Node, Vertebra

class Exon:
    def __init__(self, transcript, beg, end, seq):
        self.chrom = transcript.chrom
        self.strand = transcript.strand
        self.beg = int(beg)
        self.end = int(end)
        self.seq = seq
        self._transcript = weakref.proxy(transcript)
        self._vertebras = self._makeBackbone()
        self._dirty_references = []
        
        if not len(seq):
            sys.exit("ERROR: no exon sequence - {}\n".format(' '.join((chrom, strand, str(beg), str(end)))))
        elif len(seq) != self.end - self.beg + 1:
            sys.exit("ERROR: wrong exon sequence length - {}\n".format(' '.join((chrom, strand, str(beg), str(end)))))
    
    def _makeBackbone(self):
        vertebras = [] # exon backbone
        samples = self._transcript.getSamples()
        for i, nucl in enumerate(self.seq):
            pos = self.beg + i
            vertebra = Vertebra(pos, nucl, samples)
            if len(vertebras):
                if self.strand == '+':
                    vertebras[-1].setNext(vertebra)
                else:
                    vertebra.setNext(vertebras[-1])
            vertebras.append(vertebra)
        return vertebras
    
    def _snp_validation(self, snp):
        if snp.beg < self.beg and self.beg <= snp.end or snp.beg <= self.end and self.end < snp.end:
            print("WARNING: SNP intersects with exon borders - SNP {} ...skipped".format(snp.id))
            return (None, None)
        if snp.end < self.beg or self.end < snp.beg:
            sys.exit("ERROR: SNP position is out of exon - SNP {} ...skipped".format(snp.id))
        
        if self.strand != snp.strand:
            snp.setComplement()
        snp_start = snp.beg - self.beg
        snp_stop = snp.end - self.beg
        if snp_start < 0 or len(self._vertebras) < snp_stop + 1:
            sys.exit("ERROR: exon modification is out of the exon boundary - {}\n".format(snp.id))
        
        snp_ref = self.seq[snp_start : snp_stop + 1]
        if snp_ref.upper() != snp.ref.seq.upper():
            print("WARNING: wrong reference allele ({}<>{}) - SNP {} ...skipped".format(snp.ref.seq, snp_ref, snp.id))
            return (None, None)
        
        if self.strand == '-':
            (snp_start, snp_stop) = (snp_stop, snp_start)
        return (snp_start, snp_stop)
    
    def getFirstVertebra(self):
        if len(self._vertebras):
            return self._vertebras[0]
        return None
    
    def getLastVertebra(self):
        if len(self._vertebras):
            return self._vertebras[-1]
        return None

    def modify(self, snp):
        vertebras = self._vertebras

        (snp_exbeg, snp_exend) = self._snp_validation(snp)
        if snp_exbeg is not None and snp_exend is not None:
            self._transcript.appendVariation(snp, vertebras[snp_exbeg], vertebras[snp_exend]) # save variation backbone to find nonsynonymous SNP
            no_reference = True
            for allele in snp.getSampleAlleles():
                allele_seq = list(reversed(allele.seq)) if self.strand == '-' else allele.seq
                if allele.isReference(): # path for the reference allele
                    for allele_pos, pos in enumerate(range(snp_exbeg, snp_exend - 1, -1) if self.strand == '-' else range(snp_exbeg, snp_exend + 1)):
                        ref_node = vertebras[pos].getRefNode()
                        ref_node.appendAllele(allele)
                        if allele_seq[allele_pos] != ref_node.nucl:
                            sys.exit("ERROR: wrong exon modification - reference allele mismatch for SNP {}\n".format(snp.id))
                        vertebras[pos].addNodeToGraph(ref_node)
                    for sample in self._transcript.getSamples():
                        if not allele.getSample(sample):
                            self._dirty_references.append({'pos': snp_exbeg, 'sample': sample})
                    no_reference = False
                else: # path for alternative alleles
                    start_node = None
                    end_node = None
                    
                    for pos, nucl in enumerate(allele_seq):
                        pos = snp.beg + pos
                        pos_overhead = 0
                        if snp.end < pos:
                            pos_overhead = pos - snp.end
                            pos = snp.end
                        new_node = Node(pos, pos_overhead, nucl, allele)
                        if end_node:
                            end_node.setNext([new_node])
                        else:
                            start_node = new_node
                        end_node = new_node
                    
                    vertebras[snp_exbeg].addNodeToGraph(start_node) # join the allele start to the graph
                    vertebra = vertebras[snp_exend].getNext()
                    if vertebra:
                        end_node.setNext(vertebra.getNodes()) # join the allele end to the graph
            if no_reference:
                for sample in self._transcript.getSamples():
                    self._dirty_references.append({'pos': snp_exbeg, 'sample': sample})
        return None
    
    def cleanup(self): # remove harmful reference paths (arise with overlapped variations)
        for dirty_ref in self._dirty_references:
            vertebra = self._vertebras[dirty_ref['pos']]
            ref_node = vertebra.getRefNode()
            ref_node.removeSample(dirty_ref['sample'])
# end of Exon
