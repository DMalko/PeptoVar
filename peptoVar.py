#!/usr/bin/env python3

import sys
import re
import argparse
import time
import datetime

try:
    from pysam import VariantFile, VariantHeader
except:
    sys.exit("Can't parse VCF file.\nPlease, install pysam module for python3\n")
from lib.gff import Gff
from lib.transcript import Transcript
from lib.output import OutFileContainer, ColorPrint
from lib.input import InFile
from lib.variation import SNP
from lib.seqtools import PeptComparator
from lib.comparator import UniPep

# -debug -peptlen 8 -var nonsyn -gff ./data/Homo_sapiens.GRCh38.85.chromosome.10.gff3 -fasta ./data/Homo_sapiens.GRCh38.dna.chromosome.10.fa  -vcf ./data/ALL.chr10.vcf.gz
# -debug -peptlen 0 -var nonsyn -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz
# -trnlist transcript:ENST00000368799 -samples HG00242 NA20588 -peptlen 9 -var nonsyn -gff ./ensembl1/chr1.gff3 -vcf ./ensembl1/chr1.vcf.gz -fasta ./ensembl1/chr1.fa
# -samples T00001 T00002 -peptlen 9 -var nonsyn -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz

DEBUG = 0
if DEBUG:
    from lib.seqmod import Sequence

def main():
    start_time = time.time()
    total_tr = 0
    cprint = ColorPrint()
    
    input_parser = argparse.ArgumentParser(description='PeptoVar - Peptides on Variations: program for personalization of protein coding genes and peptidomes generation.')
    input_parser.add_argument('-gff', metavar='file.gff', default=None, help='GFF input file', required=False)
    input_parser.add_argument('-fasta', metavar='file.fasta', default=None, help='FASTA input file', required=False)
    input_parser.add_argument('-seq', metavar='file.data', default=None, help='DATA input file', required=False)
    input_parser.add_argument('-vcf', metavar='file.vcf.gz', default=None, help='bgzip-compressed VCF input file (need index file)', required=False)
    input_parser.add_argument('-tmpdir', metavar='dirpath', default=None, help='TEMP directory', required=False)
    input_parser.add_argument('-samples', metavar='name', nargs='+', default=list(), help='sample name or pair of names in VCF file; for two samples (donor/recipient) only unique peptides will be represented)', required=False)
    input_parser.add_argument('-minaf', metavar='THRESHOLD', type=float, default=0, help='allele frequency (AF) threshold; alleles with AF < THRESHOLD will be ignored (AF=0 will be set for alleles with no data)', required=False)
    input_parser.add_argument('-var', metavar='all | nonsyn', help='save translated polymorphisms (all or only non-synonymous)')
    input_parser.add_argument('-nopt', action='store_false', default=True, help='do not use optimization (may cause high CPU load and memory usage)')
    input_parser.add_argument('-peptlen', metavar='LENGTH', nargs='+', type=int, default=list(), help='lengths of peptides (0 - full-length proteins)', required=False)
    input_parser.add_argument('-outdir', metavar='dirpath', default='./output', help='output directory (will be created if not exists, default=./output)', required=False)
    input_parser.add_argument('-indir', metavar='dirpath', default=None, help='input directory for files *.vcf, *.gff and *.fasta - if no sequence in GFF file; the files MUST have the same name for each locus (chromosome)', required=False)
    input_parser.add_argument('-trnlist', metavar='transcriptID', nargs='+', default=list(), help='list of transcriptID for processing', required=False)
    input_parser.add_argument('-trnfile', metavar='transcriptID.txt', default=None, help='one column text file with transcriptID for processing', required=False)
    
    args = input_parser.parse_args()
    gff_file = args.gff
    fasta_file = args.fasta
    seq_file = args.seq
    vcf_file = args.vcf
    tmp_dir = args.tmpdir
    save_var = args.var
    optimization = args.nopt
    min_af = args.minaf
    trnlist = args.trnlist
    trnfile = args.trnfile
    pept_len = args.peptlen
    pept_len.sort()
    do_prot = False
    if len(pept_len) and pept_len[0] == 0:
        do_prot = True
        pept_len.pop(0)
    outdir = args.outdir
    indir = args.indir
    samples = args.samples
    
    peptdb = None
    protdb = None
    outfiles = OutFileContainer(outdir)
    
    if not(vcf_file or indir):
        cprint.printFail("\nYou have no input data!\n")
        input_parser.print_help()
        exit()
    
    transcript_set = {}
    if trnfile:
        with open(trnfile, "r") as trnlistfile:
            for line in trnlistfile:
                trnid = line.strip()
                if len(trnid) > 0:
                    transcript_set[trnid] = 0
        trnlistfile.close()
    for trnid in trnlist:
        transcript_set[trnid] = 0
    
    if len(samples) == 0:
        samples.append('virtual')
        cprint.printWarning('MODE: virtual sample')
        peptdb = UniPep(samples[0], '-', tmp_dir)
        protdb = UniPep(samples[0], '-', tmp_dir)
    elif len(samples) == 1:
        cprint.printWarning('MODE: sample {}'.format(samples[0]))
        peptdb = UniPep(samples[0], '-', tmp_dir)
        protdb = UniPep(samples[0], '-', tmp_dir)
    elif len(samples) == 2:
        cprint.printWarning('MODE: transplantation')
        peptdb = UniPep(samples[0], samples[1], tmp_dir)
        protdb = UniPep(samples[0], samples[1], tmp_dir)
    else:
        cprint.printFail("\nCan't take more than two samples\n")
        exit()
    
    input_bulk = []
    if vcf_file:
        if gff_file:
            input_bulk.append({'name': 'input', 'vcf': vcf_file, 'gff': gff_file, 'fasta': fasta_file, 'seq': seq_file})
        else:
            vcf_file = re.sub('.*\/', '', vcf_file)
            cprint.printWarning("No GFF file for VCF {} ...skipped!".format(vcf_file))
    if indir:
        infiles = InFile(indir)
        input_bulk.extend(infiles.bulk)
    
    if DEBUG:
        errlog = open("ERROR.log", 'w')
    
    for fileset in input_bulk:
        print("{} files parsing...\n".format(fileset['name']))
        vcf = VariantFile(fileset['vcf'])
        gff = Gff(fileset['gff'], tmp_dir)
        if fileset['seq']:
            gff.attachSeq(fileset['seq'])
        else:
            gff.attachSeq(fileset['fasta'])
        
        for trn_id in gff.getTranscriptsID():
            if len(transcript_set):
                if trn_id not in transcript_set:
                    continue
                transcript_set[trn_id] = 1
            
            total_tr += 1
            transcript = Transcript(trn_id)
            transcript.setSamples(samples)
            
            debugseq = None
            for exon in gff.getTranscriptExons(trn_id):
                transcript.appendExon(exon.chrom, exon.strand, exon.beg, exon.end, exon.seq)
                if DEBUG:
                    if not debugseq:
                        debugseq = Sequence(transcript.chrom, transcript.id, transcript.getSamples())
                    debugseq.append(exon.strand, exon.beg, exon.end, exon.seq)

            cprint.printOk("IN PROCESS #{}: {} (locus {})".format(total_tr, transcript.id, transcript.chrom))
            
            n_snp = 0
            print("transript modification...") #print('\x1b[2K\r')
            for exon in transcript.getExons():
                for var in vcf.fetch(exon.chrom, exon.beg - 1, exon.end):
                    try:
                        snp = SNP(var, samples)
                    except:
                        cprint.printFail("\nWrong sample(s) name\n")
                        exit()
                    else:
                        if snp.filterAF(min_af):
                            if DEBUG:
                                debugseq.modify(snp)
                            exon.modify(snp)
                            n_snp += 1
                exon.cleanup()
            
            print("graph building...")
            transcript.makeGraph()
            trnVariations = transcript.translateVariations()
            
            if optimization:
                print("graph optimization...")
                if save_var == 'nonsyn':
                    trnVariations = transcript.joinSynonymPathes()
                else:
                    transcript.joinSynonymPathes()
            
            if save_var: # save translated alleles in file
                print("variation processing...")
                for trn_var in trnVariations:
                    outfiles.writeVariation(transcript.id, trn_var)
            
            if do_prot:
                print("protein processing...")
                sample_proteins = transcript.getProteins(optimization)
                for sample_name in sample_proteins:
                    protdb.load(sample_name, sample_proteins[sample_name])
            
            if len(pept_len):
                print("peptide processing...")
                sample_peptides = transcript.getPeptides(pept_len, optimization)
                for sample_name in sample_peptides:
                    peptdb.load(sample_name, sample_peptides[sample_name]) # 'chrom', 'transcript_id', 'sample', 'allele1', 'allele2', 'beg', 'end', 'fshifts_before', 'variations(positions_in_matrix)', 'peptide', 'matrix'
                    if DEBUG:
                        debugseq.add2compare(sample_peptides[sample_name])
            
            if DEBUG:
                debugseq.translate()
                debugseq.peptides(pept_len)
                for misspept in debugseq.compare():
                    errlog.write(misspept)
            
            print("{} ...ok\n".format(transcript.id))
    
    print("writing to output...")
    if len(pept_len):
        peptdb.flush()
        peptdb.doIndex()
        for sample in samples:
            if len(samples) > 1:
                output = peptdb.getUnique(sample)
            else:
                output = peptdb.getAll(sample)
            for rec in output:
                outfiles.writePeptide(rec)
        peptdb.close()
    
    if do_prot:
        protdb.flush()
        protdb.doIndex()
        for sample in samples:
            if len(samples) > 1:
                output = protdb.getUnique(sample)
            else:
                output = protdb.getAll(sample)
            for rec in output:
                outfiles.writeProtein(rec)
        protdb.close()
    
    notfound = set()
    for trnid, count in transcript_set.items():
        if not count:
            notfound.add(trnid)
    if len(notfound):
        for trnid in notfound:
            outfiles.writeWarning(trnid)
        cprint.printWarning("\nWARNING: {} transcript(s) not found (see warnings.txt file)\n".format(len(notfound)))
    
    stop_time = time.time()
    cprint.printOk("...the job is done (execution time: {})".format(time.strftime("%H:%M:%S", time.gmtime(stop_time - start_time))))
# end of main

if __name__ == '__main__':
    main()
