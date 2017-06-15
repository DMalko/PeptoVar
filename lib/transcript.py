#!/usr/bin/env python3

import sys
from lib.exon import Exon
from lib.sample import Sample
from lib.graph import Node, FShiftPath, CodonEmptyPrefix, CodonPrefix, CodonEmptySuffix
from lib.trnvariation import TrnAllele, TrnVariation
from lib.peptide import Peptide, PeptideDriver
from lib.seqtools import translate, complement, is_stopcodon

class Transcript:
    def __init__(self, mrna_id , chrom = None, strand = None):
        self.id = mrna_id
        self.chrom = chrom
        self.strand = strand
        self._exon_number = 0
        self._exons = [] # to handle variations at the end of CDS the last exon must include the STOP CODON !!!!
        self._backbone_head = None # the head vertebra of the backbone
        self._graph_start = None # start nodes of the graph
        self._used_nodes = [] # array of used graph nodes
        self._variations = [] # array of variation linked to the backbones
        self._samples = [] # all samples in the run
        self._codon_start_nodes = [] # array of the first nodes in codons (for each sample separately)
    
    def cleanUsed(self):
        for used_node in self._used_nodes:
            used_node.cleanUsed() # graph clean up
        nodes = self._used_nodes
        self._used_nodes = []
        return nodes
    
    def setSamples(self, sample_names):
        for name in sample_names:
            if name == 'virtual':
                self._samples.append(Sample(name, 1, 1))
            else:
                self._samples.append(Sample(name, 1, 0))
                self._samples.append(Sample(name, 0, 1))
    
    def getSamples(self):
        return self._samples
    
    def appendExon(self, chrom, strand, beg, end, seq):
        beg = int(beg)
        end = int(end)
        if not self.chrom:
            self.chrom = chrom
        if not self.strand:
            self.strand = strand
        
        if beg < 1 or end < 1:
            sys.exit("ERROR: mRNA failure (wrong exon location) - {}\n".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        if end < beg:
            sys.exit("ERROR: mRNA failure (wrong exon orientation) - {}\n".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        if len(seq) != end - beg + 1:
            sys.exit("ERROR: mRNA failure (wrong exon sequence length) - {}\n".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        if len(self._exons):
            if beg <= self._exons[-1].end:
                sys.exit("ERROR: mRNA failure (wrong exon priority) - {}\n".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        
        new_exon = Exon(self, beg, end, complement(seq) if strand == '-' else seq)
        if len(self._exons):
            # tie the new exon with other to get the continuous matrix backbone
            if self.strand == '+':
                vert = self._exons[-1].getLastVertebra()
                vert.setNext(new_exon.getFirstVertebra())
            else:
                vert = new_exon.getFirstVertebra()
                vert.setNext(self._exons[-1].getLastVertebra())
                self._backbone_head = new_exon.getLastVertebra()
        else:
            if self.strand == '+':
                self._backbone_head = new_exon.getFirstVertebra()
            else:
                self._backbone_head = new_exon.getLastVertebra()
        self._exons.append(new_exon)
        
        # check data consistency #
        vert = new_exon.getLastVertebra() if self.strand == '-' else  self._exons[0].getFirstVertebra()
        vert_num = 0
        while vert:
            vert_num += 1
            vert = vert.getNext()
        mtx_len = 0
        for exon in self._exons:
            mtx_len += exon.end - exon.beg + 1
        if vert_num != mtx_len:
            sys.exit("ERROR: faulty node sequence - {}\n".format(self.id))
        
        self._exon_number += 1
        return self._exon_number
    
    def getExons(self):
        return self._exons
    
    def getVariations(self):
        return self._variations
    
    def appendVariation(self, snp, beg_vertebra, end_vertebra):
        self._variations.append(TrnVariation(snp, beg_vertebra, end_vertebra))
        return len(self._variations)
    
    def makeGraph(self):
        if not self._backbone_head:
            return None
        elif self._graph_start:
            return self._graph_start
        
        vertebra = self._backbone_head
        self._graph_start = Node(0, 0, '')
        self._graph_start.setNext(vertebra.getNodes())
        while vertebra:
            vertebra.addNodeToGraph()
            vertebra = vertebra.getNext()
        
        for sample in self._samples:
            self._findFrameShifts(sample)
            nodes = self._attachPrefixes(sample)
            self._codon_start_nodes.append((sample, nodes))
        return self._codon_start_nodes
    
    def _findFrameShifts(self, sample):
        def pathFinder(path):
            node = path['node']
            frame = path['frame']
            fshift_path = path['fsh_path']
            next_path = []
            if len(node.getNext(sample)):
                frame += 1
                if frame > 2: frame = 0
                
                if node.isFrameShift():
                    allele_id = node.getAlleleID()
                    fshift_path = fshift_path.clonePath()
                    fshift_path.appendAlleleID(allele_id)
                for next_node in node.getNext(sample):
                    if frame == 0:
                        if next_node.appendFShiftPath(fshift_path): # check if the path has not been used yet (STOP if it's used)
                            next_path.append({'node': next_node, 'frame': frame, 'fsh_path': fshift_path})
                    else:
                        next_path.append({'node': next_node, 'frame': frame, 'fsh_path': fshift_path})
            return next_path
        # end of pathFinder()
        
        tree = []
        fshift_startpath = FShiftPath(sample)
        for node in self._graph_start.getNext(sample):
            node.appendFShiftPath(fshift_startpath)
            tree.append({'node': node, 'frame': 0, 'fsh_path': fshift_startpath})
        while len(tree):
            new_tree = []
            for path in tree:
                new_path = pathFinder(path)
                new_tree.extend(new_path)
            tree = new_tree
    
    def _attachPrefixes(self, sample):
        def prefixDriver(node1): # node1 is the first nucleotide in the codon
            next_nodes = []
            if node1.isUsed():
                return next_nodes # STOP - the node have been used with other path
            if not node1.getSample(sample):
                return next_nodes # the node is not in the sample path
            for node2 in node1.getNext(sample): # the second nucleotide
                for node3 in node2.getNext(sample): # the third nucleotide
                    node1.attachPrefix(CodonEmptyPrefix(sample, node1))
                    node2.attachPrefix(CodonPrefix(sample, [node1]))
                    node3.attachPrefix(CodonPrefix(sample, [node1, node2]))
                    for vertebra in (node1.getVertebra(), node2.getVertebra(), node3.getVertebra()):
                        if vertebra:
                            vertebra.appendSample(sample)
                    if is_stopcodon(node1.nucl + node2.nucl + node3.nucl):
                        continue # there is no reason to attach prefixes to untranslated codons
                    for next_node in node3.getNext(sample):
                        next_nodes.append(next_node)
            node1.setUsed() # mark the node as `used`
            self._used_nodes.append(node1) # save used nodes to clean the marks before the next run
            
            return next_nodes
        
        tree = self._graph_start.getNext(sample)
        while len(tree):
            new_tree = []
            for node in tree:
                new_nodes = prefixDriver(node)
                new_tree.extend(new_nodes)
            tree = new_tree
        return self.cleanUsed() # array of the first nodes in codons can be used to make peptides 
    
    def translateVariations(self):
        for var in self._variations:
            beg_vertebra = var.beg_vertebra
            end_vertebra = var.end_vertebra.getNext()
            for allele in var.snp.getSampleAlleles():
                allele_seq = allele.seq if self.strand == '+' else allele.seq[::-1]
                for sample in self._samples:
                    if beg_vertebra.getSample(sample): # take only translated alleles
                        prefixes = beg_vertebra.getPrefixes(sample)
                        if end_vertebra:
                            siffixes = end_vertebra.getSuffixes(sample) # the last exon must have STOP CODON to get suffix from the last codon
                        else:
                            siffixes = [CodonEmptySuffix(sample)]
                        for prefix in prefixes: # get allele translations for each pair prefix/suffix
                            used_suffix_seq = set()
                            for suffix in siffixes:
                                suffix_len = 3 - (len(prefix.seq) + len(allele.seq)) % 3
                                if suffix_len == 3: suffix_len = 0
                                suffix = suffix.getTrimmed(suffix_len)
                                if suffix.seq in used_suffix_seq:
                                    continue
                                used_suffix_seq.add(suffix.seq)
                                mtx = prefix.seq + allele_seq + suffix.seq
                                if len(mtx) < 3:
                                    continue
                                (trn, tail) = translate(mtx)
                                var.appendTrnAllele(TrnAllele(allele, prefix, suffix, mtx, trn))
        return self._variations
    
    def getProteins(self, optimization = False):
        return self.getPeptides([0], optimization)
    
    def getPeptides(self, pept_len, optimization = False):
        pept_container = {}
        codon_start_nodes = []
        if pept_len and len(pept_len) > 0:
            if len(pept_len) == 1 and pept_len[0] == 0: # proteins
                for sample in self._samples:
                    codon_start_nodes.append((sample, self._graph_start.getNext(sample, '')))
            elif self._graph_start: # peptides
                codon_start_nodes = self._codon_start_nodes
        
        for sample, start_nodes in codon_start_nodes:
            for node in start_nodes:
                peptDriver = PeptideDriver(self.chrom, self.id, sample, pept_len, optimization)
                peptDriver.appendNode(node)
                peptDriver.setPeptContainer(pept_container)
                driverArray = [peptDriver]
                while driverArray:
                    new_driverArray = []
                    for driver in driverArray:
                        nodes = driver.getNext()
                        if len(nodes) > 1:
                            for node in nodes:
                                new_peptDriver = driver.copy() # !!!
                                if new_peptDriver.appendNode(node):
                                    new_driverArray.append(new_peptDriver)
                        elif len(nodes) == 1:
                            if driver.appendNode(nodes[0]):
                                new_driverArray.append(driver)
                    driverArray = new_driverArray
        return pept_container
    
    def getBackboneProtein(self):
        mtx = ""
        vertebra = self._backbone_head
        while vertebra:
            mtx += vertebra.getRefNode().nucl
            vertebra = vertebra.getNext()
        return translate(mtx)
    
    def joinSynonymPathes(self):
        nonsynonym_var = []
        for trn_var in self._variations:
            trn_var.joinSynonymAlleles()
            if trn_var.isNonSynonymous():
                nonsynonym_var.append(trn_var)
        return nonsynonym_var
# end of Transcript
