#!/usr/bin/env python3

import errno
import os
import re
import tempfile

class ColorPrint:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
    def printWarning(self, text):
        print(self.WARNING + text + self.ENDC)
    
    def printOk(self, text):
        print(self.OKGREEN + text + self.ENDC)
        
    def printFail(self, text):
        print(self.FAIL + text + self.ENDC)

class OutFile:
    def __init__(self, dir_name, file_name):
        if re.search(r"\/$", dir_name) is None:
            dir_name += '/'
        self._file_path = dir_name + file_name
        self._file_handle = open(self._file_path, "w")
        self._written = 0
    
    def __del__(self):
        if self._file_handle and not self._file_handle.closed:
            self._file_handle.close()
    
    def writeCSV(self, data_array):
        if self._written:
            self._file_handle.write("\n")
        else:
            self._written = 1
        self._file_handle.write("\t".join(str(item) for item in data_array))
        
    def close(self):
        self.__del__()
# end of OutFile

class OutFileContainer:
    def __init__(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
        self._path = path
        self._files_prot = {}
        self._files_pept = {}
        self._file_var = None
        self._file_warn = None
    
    def writeWarning(self, warn):
        if not self._file_warn:
            self._file_warn = OutFile(self._path, 'warnings.txt')
        self._file_warn.writeCSV((warn))
    
    def writeProtein(self, protein_rec):
        sample_name = protein_rec[2]
        if sample_name not in self._files_prot:
            self._files_prot[sample_name] = OutFile(self._path, sample_name + '.prot.csv')
            self._files_prot[sample_name].writeCSV(('chrom', 'transcript_id', 'sample', 'allele1', 'allele2', 'beg', 'end', 'variations(positions_in_matrix)', 'protein', 'matrix'))
        outfile = self._files_prot[sample_name]
        outfile.writeCSV(protein_rec[0:6] + protein_rec[7:])
    
    def writePeptide(self, peptide_rec):
        sample_name = peptide_rec[2]
        if sample_name not in self._files_pept:
            self._files_pept[sample_name] = OutFile(self._path, sample_name + '.pept.csv')
            self._files_pept[sample_name].writeCSV(('chrom', 'transcript_id', 'sample', 'allele1', 'allele2', 'beg', 'end', 'fshifts_before', 'variations(positions_in_matrix)', 'peptide', 'matrix'))
        outfile = self._files_pept[sample_name]
        outfile.writeCSV(peptide_rec)
    
    def writeVariation(self, transcript_id, trn_var):
        if not self._file_var:
            self._file_var = OutFile(self._path, 'variations.csv')
            self._file_var.writeCSV(('transcript_id', 'variation_id', 'beg', 'end', 'allele_id', 'sample', 'sample_allele_1', 'sample_allele_2', 'nonsynonymous', 'fshifts_before', 'prefix', 'allele', 'suffix', 'translation'))
        outfile = self._file_var
        
        var_id = trn_var.id
        beg = trn_var.beg_vertebra.pos
        end = trn_var.end_vertebra.pos
        trn_alleles = trn_var.getTrnAlleles()
        
        for trn_allele_id in trn_alleles:
            nonsyn = 0
            if trn_var.isNonSynonymous(trn_allele_id):
                nonsyn = 1
            for trn_allele in trn_alleles[trn_allele_id]:
                prefix = trn_allele.prefix
                suffix = trn_allele.suffix
                if prefix.sample == suffix.sample and trn_allele.allele.getSample(prefix.sample): # a sample allele must be in the native context!!!
                    prefix_fsh = []
                    for fshift in prefix.getFShiftPathSet():
                        prefix_fsh.append(fshift.id)
                    
                    outfile.writeCSV((transcript_id, var_id, str(beg), str(end), trn_allele.id, prefix.sample.name, str(prefix.sample.allele_1), str(prefix.sample.allele_2), str(nonsyn), ';'.join(prefix_fsh), prefix.seq, trn_allele.allele.seq, suffix.seq, trn_allele.trn))
    
    def close(self):
        for file in list(self._files_prot.values()):
            file.close()
        for file in list(self._files_pept.values()):
            file.close()
        if self._file_var:
            self._file_var.close()
        if self._file_warn:
            self._file_warn.close()
# end of OutFiles


    