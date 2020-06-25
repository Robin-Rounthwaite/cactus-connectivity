#%%

"""
I want toilified, parallelized pipeline that calculates all_to_all_liftovers.

Note: for now, I'm manually crafting my "assembly_files" dict from the seq.txt used in cactus. Shouldn't be too hard to automate that, though. 
        Maybe even extract the info from the cactus graph itself, if I"m feeling ambitious (and it records the original fasta files it was made from).
"""

from toil.common import Toil
from toil.job import Job
import os
import subprocess
from Bio import SeqIO

import collections as col

#first step is to make the full_beds.
def get_contig_lengths(assembly):
    lengths = dict()
    asm = SeqIO.index(assembly, "fasta")
    for contig_id, seq in asm.items():
        lengths[contig_id] = len(seq)
    return lengths

def write_full_bed(contig_lengths, out_bed):
    with open(out_bed, "w") as outf:
        for contig_id, length in contig_lengths.items():
            outf.write(contig_id + "\t" + "0" + "\t" + str(length) + "\n")
    return out_bed

#Second step is to call liftover on each possible combination of assembly.
def liftover(hal_file, source_assembly, source_full_bed, target_assembly, out_bed):
    #todo: edit:
    subprocess.call(["halLiftover",  hal_file, source_assembly, source_full_bed, target_assembly, out_bed])
    
def all_to_all_liftovers(assembly_files, hal_file, output_dir):
    """assembly_files is a dict with key: assembly name and value: assembly_file.
    """
    # first, calculate lengths of contigs in each asm:
    lengths = dict()
    for asm, asm_file in assembly_files.items():
        lengths[asm] = get_contig_lengths(asm_file)
    
    # Then, make the full.bed files, which will act as srcBed in the liftover. This way,
    # the liftover will look for where the target genome is mapped to all possible locations
    # in the src genome.
    full_beds = dict()
    for asm in assembly_files:
        #todo: change tmp dir to job.fileStore.gettmpfile (in write_full_bed):
        out_bed_dir = "tmp_outbeds/"
        full_beds[asm] = write_full_bed(lengths[asm], out_bed_dir + asm + ".full.bed")

    #todo: make it run all liftovers. Currently running minimal test with one liftover.
    for asm1 in assembly_files:
        for asm2 in assembly_files:
            if asm1 == asm2:
                continue 
            out_bed = output_dir +asm1 + "_source_" + asm2 + "_target_liftover.bed"
            liftover(hal_file, asm1, full_beds[asm1], asm2,  out_bed)

    
assembly_dir = "./asms/"
output_dir = "./liftovers/"
assembly_files = {"HG03098_paf_chr21": assembly_dir + "HG03098_paf_chr21.fa", "HG03492_paf_chr21": assembly_dir + "HG03492_paf_chr21.fa", "hg38_chr21": assembly_dir + "hg38_chr21.fa"}
hal_file = "ref_based_small_chr21.hal"
all_to_all_liftovers(assembly_files, hal_file, output_dir)