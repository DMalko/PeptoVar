#!/usr/bin/env python3

import re
import os
import random
import string
import tempfile

class Seq:
    def __init__(self, seq_id, tmp_dir = None):
        self.id = seq_id
        self._file = tempfile.TemporaryFile(mode='w+', suffix='.seq', dir=tmp_dir)
        self.len = 0
        self._shiftpos = 0
    
    def __del__(self):
        if self._file and not self._file.closed:
            self._file.close()
    
    def get(self, beg, end):
        if self.len < end or end < beg:
            return None
        self._file.seek(self._shiftpos + beg - 1)
        return self._file.read(end - beg + 1)
    
    def append(self, seq):
        self._file.write(seq)
        self.len += len(seq)
        
    def set(self, fh, pos = 0):
        self._file = fh
        self._shiftpos = pos
        fh.seek(0, 2)
        self.len = fh.tell() - pos
        fh.seek(0)
        
    def flush(self):
        self._file.flush()
# end of Seq

class Fasta:
    def __init__(self, tmp_dir = None):
        self._seq = {}
        self._tmp = tmp_dir
    
    def appendSeq(self, file):
        file = open(file, 'r')
        pos = 0
        seq_id = ''
        for line in file:
            pos = len(line)
            re_header = re.match(">(\S+)", line)
            if re_header:
                seq_id = re_header.group(1)
            else:
                print("ERROR: wrong format of sequence file")
                exit()
            break
        seq = Seq(seq_id, self._tmp)
        seq.set(file, pos)
        self._seq[seq_id] = seq
    
    def appendFasta(self, file, pos = 0):
        file = open(file, 'r')
        file.seek(pos)
        seq = None
        hdr_pattern = re.compile('>(\S+)')
        n = 0
        for line in file:
            n += 1
            if line[0] == '>':
                match = re.match(hdr_pattern, line)
                if match:
                    seq_id = match.group(1)
                    if seq: seq.flush()
                    seq = Seq(seq_id, self._tmp)
                    self._seq[seq_id] = seq
                else:
                    sys.exit("ERROR: wrong FASTA file format (line in FASTA: {} )\n".format(n))
            elif seq:
                seq.append(line.strip())
        if seq:
            seq.flush()
        file.close()
    
    def getSeq(self, seq_id, beg, end):
        if seq_id in self._seq:
            return self._seq[seq_id].get(beg, end)
        return None
# end of Fasta
