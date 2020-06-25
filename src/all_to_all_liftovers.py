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
from argparse import ArgumentParser
import collections as col

def empty(job):
    """
    An empty job, for easier toil job organization.
    """
    return

#first step is to make the full_beds.
def get_contig_lengths(job, assembly):
    lengths = dict()
    asm = SeqIO.index(job.fileStore.readGlobalFile(assembly), "fasta")
    for contig_id, seq in asm.items():
        lengths[contig_id] = len(seq)
    return lengths

def write_full_bed(job, contig_lengths):
    out_bed = job.fileStore.getLocalTempFile()
    
    with open(out_bed, "w") as outf:
        for contig_id, length in contig_lengths.items():
            outf.write(contig_id + "\t" + "0" + "\t" + str(length) + "\n")

    return job.fileStore.writeGlobalFile(out_bed)

#Second step is to call liftover on each possible combination of assembly.
def liftover(job, hal_file, source_assembly, source_full_bed, target_assembly):
    out_bed_tmp = job.fileStore.getLocalTempFile()
    out_bed = job.fileStore.writeGlobalFile(out_bed_tmp)

    subprocess.call(["halLiftover", job.fileStore.readGlobalFile(hal_file), source_assembly, job.fileStore.readGlobalFile(source_full_bed), target_assembly, job.fileStore.readGlobalFile(out_bed)])
    
    return out_bed
    
def all_to_all_liftovers(job, assembly_files, assembly_lengths, hal_file):
    """assembly_files is a dict with key: assembly name and value: assembly_file.
    """
    lengths_jobs = job.addFollowOnJobFn(empty)
    
    # Then, make the full.bed files, which will act as srcBed in the liftover. This way,
    # the liftover will look for where the target genome is mapped to all possible locations
    # in the src genome.
    full_beds = dict()
    for asm in assembly_files:
        full_beds[asm] = lengths_jobs.addChildJobFn(write_full_bed, assembly_lengths[asm]).rv()

    full_beds_jobs = lengths_jobs.addFollowOnJobFn(empty)

    #liftovers is nested dict, with key:(target_asm), value:<dict, with key:source_asm, value:<list of liftover_files with target_asm as target> >
    liftovers = dict()
    for target_asm in assembly_files:
        liftovers[target_asm] = dict()
        for source_asm in assembly_files:
            
            if source_asm == target_asm:
                continue

            # liftovers[target_asm][source_asm] = full_beds_jobs.addChildJobFn(liftover, hal_file, source_asm, full_beds[source_asm], target_asm, cores=1).rv()
            liftovers[target_asm][source_asm] = full_beds_jobs.addChildJobFn(liftover, hal_file, source_asm, full_beds[source_asm], target_asm).rv()
    
    return liftovers

def main():
    # if I wanted to make this into a true command line tool, I'd fill out the parser.
    # Instead, I'm just going to add the bare minimum for making a workflow. 
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    # parser.add_argument(
    #     'output_dir', help='The dir to save the output, target bedfiles.', type=str)
    options = parser.parse_args()

    assembly_dir = "./asms/"
    output_dir = "./liftovers/"
    assembly_files = {"HG03098_paf_chr21": assembly_dir + "HG03098_paf_chr21.fa", "HG03492_paf_chr21": assembly_dir + "HG03492_paf_chr21.fa", "hg38_chr21": assembly_dir + "hg38_chr21.fa"}
    hal_file = "ref_based_small_chr21.hal"

    with Toil(options) as workflow:
        if not workflow.options.restart:
            #importing files:
            for asm, asm_file in assembly_files.items():
                assembly_files[asm] = workflow.importFile("file://" + os.path.abspath(asm_file))
            
            hal_file = workflow.importFile("file://" + os.path.abspath(hal_file))
                
            #todo: update here, for running not in cactus_connectivity, need assembly_lengths, possibly other things.
            liftovers = workflow.start(Job.wrapJobFn(all_to_all_liftovers, assembly_files, hal_file, output_dir))
            # liftovers = workflow.start(Job.wrapJobFn(all_to_all_liftovers, assembly_files, hal_file, output_dir, cores=3))

            for target_asm, liftovers_dict in liftovers.items():
                for source_asm, liftover_file in liftovers_dict.items():
                    workflow.exportFile(liftover_file, 'file://' + os.path.abspath(output_dir) + "/" + source_asm + "_source_" + target_asm + "_target_liftover.bed")

        else:
            output = workflow.restart()

if __name__ == "__main__":
    main()