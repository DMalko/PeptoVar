#!/usr/bin/env python3

import sys
import weakref
from lib.seqtools import translate
from lib.sample import SampleContainer

class Node:
    def __init__(self, pos, overpos, nucl, allele = None):
        self.pos = pos
        self.overpos = overpos
        self.nucl = nucl
        self._is_ref = False if allele else True # TRUE if it belongs to the reference sequence
        self._next = [] # the next nodes
        self._vertebra = None # reference to the master vertebra
        self._allele_id = allele.id if allele else 'ref'
        self._alleles = [allele] if allele else [] # alleles of variations
        self._is_fshift = True if allele and allele.isInDel() else False # True if the node belongs to a frameshift
        self._fshift_path = uContainer() # sets of alleles with frameshifts BEFORE the node
        self._is_used = False # True if the node has been used
        self._prefix_seq = set()
        self._samples = SampleContainer()
        self._in_graph = False
        if allele:
            self.appendSamples(allele.getSamples())
    
    def setInGraph(self):
        self._in_graph = True
    
    def isInGraph(self):
        return self._in_graph
    
    def setVertebra(self, vertebra):
        self._vertebra = weakref.proxy(vertebra)
    
    def getVertebra(self):
        return self._vertebra
    
    def appendSamples(self, samples):
        self._samples.appendSamples(samples)
    
    def removeSample(self, sample):
        self._samples.removeSample(sample)
    
    def appendPrefixSeq(self, prefix_seq):
        self._prefix_seq.add(prefix_seq.lower())
        return self._prefix_seq
    
    def appendAllele(self, allele):
        if not self._is_ref: # an alternative (not reference) allele has only one allele_id
            sys.exit("ERROR: wrong allele ID for node in pos: {}\n".format(self.pos))
        self._alleles.append(allele)
        self.appendSamples(allele.getSamples())
        return len(self._alleles)
    
    def getAlleleID(self):
        return self._allele_id
    
    def getAlleles(self, sample = None):
        if sample:
            alleles = []
            for allele in self._alleles:
                if allele.getSample(sample):
                    alleles.append(allele)
            return alleles
        return self._alleles
    
    def isReference(self):
        return self._is_ref
    
    def getSample(self, sample):
        return self._samples.getSample(sample)
    
    def getPrefixSeq(self, prefix_seq):
        if prefix_seq.lower() in self._prefix_seq:
            return self._prefix_seq
        return None
    
    def getNext(self, sample = None, pefix_seq = None):
        if sample:
            next_nodes = []
            for next_node in self._next:
                if next_node.getSample(sample):
                    if pefix_seq and len(self._next) > 1: # check prefix if branches exist
                        if next_node.getPrefixSeq(pefix_seq):
                            next_nodes.append(next_node)
                    else:
                        next_nodes.append(next_node)
            return next_nodes
        return self._next
    
    def setNext(self, next_nodes):
        self._next = next_nodes
        return self._next
    
    def setUsed(self):
        self._is_used = True
        
    def cleanUsed(self):
        self._is_used = False
    
    def isUsed(self):
        return self._is_used
    
    def appendFShiftPath(self, fshift_path):
        return self._fshift_path.append(fshift_path)
    
    def getFShiftPathSet(self, sample):
        return  self._fshift_path.get(sample)
    
    def isFrameShift(self):
        return self._is_fshift
    
    def attachPrefix(self, prefix):
        if self._vertebra:
            self._vertebra.appendPrefix(prefix)
            return True
        return False
    
    def getPos(self):
        return {'pos': self.pos, 'overpos': self.overpos, 'id': self._allele_id}
# end of Node

class Vertebra:
    def __init__(self, pos, nucl, samples):
        self.pos = pos
        self._ref_node = Node(pos, 0, nucl)
        self._nodes = [] # nodes in the graph
        self._next = None # the next vertebra in a backbone
        self._prefixes = uContainer() # prefixes (CodonPrefix objects) to build the first codons for variation alleles
        self._samples = SampleContainer()
        self._all_samples = samples
    
    def appendSample(self, sample):
        self._samples.appendSample(sample)
        
    def getSample(self, sample):
        return self._samples.getSample(sample)
    
    def appendPrefix(self, prefix):
        return self._prefixes.append(prefix)
    
    def getPrefixes(self, sample):
        prefixes = self._prefixes.get(sample)
        if not len(prefixes):
            sys.exit("ERROR: no prefix in position: {}\n".format(self.pos))
        return prefixes
    
    def getSuffixes(self, sample):
        suffixes = []
        nodes = self._nodes
        if len(nodes):
            for node1 in nodes:
                if not node1.getSample(sample):
                    continue
                if len(node1.getNext(sample)):
                    for node2 in node1.getNext(sample):
                        suffixes.append(CodonSuffix(sample, [node1, node2]))
                else:
                    suffixes.append(CodonSuffix(sample, [node1, Node(0, 0, 'N')]))
        else:
            suffixes.append(CodonSuffix(sample, [Node(0, 0, 'N'), Node(0, 0, 'N')]))
        return suffixes
    
    def addNodeToGraph(self, node = None): # add alleles to the graph
        if node: # add the variation to the graph
            if node.isReference() and not len(node.getNext()) and self._next:
                node.setNext(self._next.getNodes())
            node.setVertebra(self)
            if not node.isInGraph():
                node.setInGraph()
                self._nodes.append(node)
        else: # add the reference allele to the graph if no alternatives
            node = self._ref_node
            if not len(self._nodes) and not len(node.getNext()): # if reference not in the graph
                if self._next:
                    node.setNext(self._next.getNodes())
                node.setVertebra(self)
                node.appendSamples(self._all_samples)
                if not node.isInGraph():
                    node.setInGraph()
                    self._nodes.append(node)
        return None
    
    def getRefNode(self):
        return self._ref_node
    
    def getNodes(self):
        return self._nodes
    
    def setNext(self, vertebra):
        self._next = vertebra
    
    def getNext(self):
        return self._next
# end of Vertebra

class Codon:
    def __init__(self, nodes):
        self.nodes = nodes
        self.codon = nodes[0].nucl + nodes[1].nucl + nodes[2].nuc
        self.aa = translate(self.codon)[0]
# end of Codon

class CodonPrefix:
    def __init__(self, sample, nodes = []):
        self.id = "-"
        self.seq = ""
        self.sample = sample
        self._fsh_path_set = []
        self.nodes = nodes
        if len(nodes):
            self.id = "-".join(str(id(node)) for node in nodes)
            self.seq = "".join(node.nucl.lower() for node in nodes)
            self._setFShiftPathSet(self.nodes[0])
        
    def _setFShiftPathSet(self, node):
        self._fsh_path_set = node.getFShiftPathSet(self.sample)
        
    def getAllelesID(self):
        ids = set()
        for node in nodes:
            for allele in node.getAlleles():
                ids.add(allele.id)
        return ids
        
    def getFShiftPathSet(self):
        return self._fsh_path_set
# end of CodonPrefix

class CodonEmptyPrefix(CodonPrefix):
    def __init__(self, sample, node):
        super().__init__(sample)
        super()._setFShiftPathSet(node)
# end of CodonEmptyPrefix

class CodonSuffix (CodonPrefix):
    def __init__(self, sample, nodes = []):
        super().__init__(sample, nodes)
    
    def getTrimmed(self, length = 0):
        return CodonSuffix(self.sample, self.nodes[:length])
# end of CodonSuffix

class CodonEmptySuffix(CodonSuffix):
    def __init__(self, sample):
        super().__init__(sample, [Node(0, 0, 'N'), Node(0, 0, 'N')])
# end of CodonEmptyPrefix

class uContainer:
    def __init__(self):
        self._samples = {}
        self._ids = set()
    
    def append(self, item):
        container_id = item.sample.id + item.id
        if container_id not in self._ids:
            if item.sample.id not in self._samples:
                self._samples[item.sample.id] = []
            self._samples[item.sample.id].append(item)
            self._ids.add(container_id)
            return True
        return False
    
    def get(self, sample):
        if sample.id in self._samples:
            return self._samples[sample.id]
        return []
# end of uContainer

class FShiftPath:
    def __init__(self, sample):
        self.id = "-"
        self.sample = sample
        self._path = []
    
    def appendAlleleID(self, allele_id):
        if len(self._path):
            if self._path[-1] != allele_id:
                self._path.append(allele_id)
                self.id += "," + allele_id
                return True
            else:
                return False
        else:
            self._path.append(allele_id)
            self.id = allele_id
            return True
    
    def setPath(self, path_id, path):
        self.id = path_id
        self._path = path
    
    def clonePath(self):
        cloned_path = FShiftPath(self.sample)
        cloned_path.setPath(self.id, self._path)
        return cloned_path
    
    def checkConsistency(self):
        checker = set()
        for allele_id in self._path:
            if allele_id in checker:
                return False
            else:
                checker.add(allele_id)
        return True
# end of FShiftPath
