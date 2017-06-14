#!/usr/bin/env python3

import re
import glob
import os.path
from lib.output import ColorPrint

class InFile:
    def __init__(self, dir_path):
        if os.path.isdir(dir_path):
            self.bulk = self._globbing(dir_path)
        else:
            self.bulk = []
    
    def _globbing(self, path):
        input_set = []
        re.sub(path, '', r'\/$')
        for vcf_name in glob.glob(path + "/*.vcf.gz"):
            re_name = re.match(r'(.*\/(.+))\.vcf.gz$', vcf_name)
            if re_name:
                full_name = re_name.group(1)
                name = re_name.group(2)
                gff_name = None
                gff_name1 = full_name + ".gff"
                gff_name2 = full_name + ".gff3"
                fasta_name = None
                fasta_name1 = full_name + ".fasta"
                fasta_name2 = full_name + ".fa"
                if os.path.isfile(gff_name1):
                    gff_name = gff_name1
                elif os.path.isfile(gff_name2):
                    gff_name = gff_name2
                if gff_name:
                    if os.path.isfile(fasta_name1):
                        fasta_name = fasta_name1
                    elif os.path.isfile(fasta_name2):
                        fasta_name = fasta_name2
                    input_set.append({'name': name, 'vcf': vcf_name, 'gff': gff_name, 'fasta': fasta_name, 'seq': None})
                else:
                    cprint = ColorPrint()
                    cprint.printWarning("No {}.gff or {}.gff3 file ...skipped!".format(name))
        return input_set
# end of InFile