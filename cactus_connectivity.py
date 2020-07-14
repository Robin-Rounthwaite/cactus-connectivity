# Self-made libraries:
from src import all_to_all_liftovers
from src import calculate_bases_unmapped
from src import calculate_asm_mapping_depths

from argparse import ArgumentParser
import os

from toil.common import Toil
from toil.job import Job

def get_asm_mapping_depths(job, assembly_files, hal_file):
    leader = job.addChildJobFn(all_to_all_liftovers.empty)
    
    # Part 0: calculate lengths of contigs in each asm:
    contig_lengths = dict()
    for asm, asm_file in assembly_files.items():
        contig_lengths[asm] = leader.addChildJobFn(all_to_all_liftovers.get_contig_lengths, asm_file).rv()
    lengths_jobs = leader.encapsulate()

    # Part 1: perform all_to_all_liftovers:
    liftovers = lengths_jobs.addChildJobFn(all_to_all_liftovers.all_to_all_liftovers, assembly_files, contig_lengths, hal_file).rv()

    liftovers_jobs = lengths_jobs.encapsulate()

    # Part 2: calculate the bases left unmapped on each assembly:
    mapping_depths = liftovers_jobs.addChildJobFn(calculate_asm_mapping_depths.calculate_all_mapping_depths, liftovers, contig_lengths).rv()
    mapping_depths_jobs = liftovers_jobs.encapsulate()

    #todo: change mapping_depths to a formatted output file.
    output_file = mapping_depths_jobs.addChildJobFn(asm_mapping_depths_output, mapping_depths, contig_lengths).rv()
    return output_file

def asm_mapping_depths_output(job, mapping_depths, contig_lengths):
    output_file = job.fileStore.getLocalTempFile()
    # for asm, len_dict in contig_lengths.items():
    #     print(asm, len(len_dict))
    with open(output_file, "w") as outf:
        for target_asm, (asm_mapping_depths, debug_1_if, debug_2_if) in mapping_depths.items():
            # (asm_mapping_depths, debug_1_if, debug_2_if) = mapping_depths[target_asm]

            asm_predicted_length = int()
            for depth in asm_mapping_depths:
                asm_predicted_length += asm_mapping_depths[depth]
            outf.write("sum of all bases in " + target_asm + " according to mapping_depths calc:\t" + str(asm_predicted_length) + "\ttrue length:\t" + str(get_asm_length(contig_lengths[target_asm])) + "\tratio:\t" + str(asm_predicted_length/get_asm_length(contig_lengths[target_asm])) + "\n")
            outf.write("debug_1_if "  + str(debug_1_if) +  " debug_2_if "  + str(debug_2_if) + "\n")
        asm_lengths = dict()
        for asm in contig_lengths:
            asm_lengths[asm] = get_asm_length(contig_lengths[asm])

        outf.write("\nasm_mapping_depths dictionary:\n" + str(mapping_depths) + "\n\nasm_lengths dictionary:\n" + str(asm_lengths) + "\n")


        
    return job.fileStore.writeGlobalFile(output_file)
            
def get_bases_unmapped(job, assembly_files, hal_file, options):
    leader = job.addChildJobFn(all_to_all_liftovers.empty)
    
    # Part 0: calculate lengths of contigs in each asm:
    contig_lengths = dict()
    for asm, asm_file in assembly_files.items():
        contig_lengths[asm] = leader.addChildJobFn(all_to_all_liftovers.get_contig_lengths, asm_file).rv()

    lengths_jobs = leader.encapsulate()

    # Part 1: perform all_to_all_liftovers:
    liftovers = lengths_jobs.addChildJobFn(all_to_all_liftovers.all_to_all_liftovers, assembly_files, contig_lengths, hal_file).rv()

    liftovers_jobs = lengths_jobs.encapsulate()

    # Part 2: calculate the bases left unmapped on each assembly:
    bases_unmapped = liftovers_jobs.addChildJobFn(calculate_bases_unmapped.calculate_all_bases_unmapped, liftovers, contig_lengths, options.minimum_size_gap).rv()
    bases_unmapped_jobs = liftovers_jobs.encapsulate()

    output_file = bases_unmapped_jobs.addChildJobFn(bases_unmapped_output, bases_unmapped, contig_lengths).rv()

    return output_file
    # return bases_unmapped, contig_lengths

def bases_unmapped_output(job, bases_unmapped, contig_lengths):
    output_list = list()
    output_list.append("asssembly\t(unmapped_sequence/length_of_assembly)\tunmapped_sequence\tlength_of_assembly")

    for asm, seq_unmapped in bases_unmapped.items():
        output_list.append(asm + "\t" + str(seq_unmapped/get_asm_length(contig_lengths[asm])) + "\t" + str(seq_unmapped) + "\t" + str(get_asm_length(contig_lengths[asm])))
    
    output_formatted = '\n'.join(output_list)
    
    output = job.fileStore.getLocalTempFile()
    with open(output, "w") as inf:
        inf.write(output_formatted)
    
    return job.fileStore.writeGlobalFile(output)

    
def parse_seq_file(seq_file):
    assembly_files = dict()
    with open(seq_file) as inf:
        for line in inf:
            parsed = line.split()
            assembly_files[parsed[0]] = parsed[1]
    return assembly_files

def get_asm_length(contig_lengths):
    """
    Given dict of key: contig_id, value: length of contig, returns sum of all lengths.
    """
    len_sum = int()
    for length in contig_lengths.values():
        len_sum += length
    return len_sum

def main():
    """
    Example call for small_chr21 example:
    python cactus_connectivity.py js small_chr21.txt ./halLiftover_all_to_all/ref_based_small_chr21.hal
    """
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument(
        'seq_file', help='A tab separated file with two columns: the name of each assembly, and the name of its respective fasta files in their file locations. Similar format to seqFile used as input for cactus-prepare.', type=str)
    parser.add_argument(
        'hal_file', help='The location of the hal file to be profiled.', type=str)
    #todo: minimum size gap includes more seq not mapped, or more seq as mapped?
    parser.add_argument(
        '--minimum_size_gap', help="When calculating the amount of sequence that isn't mapped, you can specify to include bases that are (not?) mapped that are smaller than the specified gapsize.TODOFINISH.", type=str)
    parser.add_argument(
        '--output', help='The dir to save the output, target bedfiles.', default='./cactus_connectivity_output.txt', type=str)
    options = parser.parse_args()
    options.minimum_size_gap = 0

    assembly_files = parse_seq_file(options.seq_file)
    # print(assembly_files)
    
    
    # # data for testing with small_chr21:
    # assembly_dir = "./halLiftover_all_to_all/asms/"
    # output = "./output/"
    # #todo: when making assembly_files+assembly_dir a command line argument, you can add it in the format of the .txt tree you send to cactus, if you like.
    # assembly_files = {"HG03098_paf_chr21": assembly_dir + "HG03098_paf_chr21.fa", "HG03492_paf_chr21": assembly_dir + "HG03492_paf_chr21.fa", "hg38_chr21": assembly_dir + "hg38_chr21.fa"}
    # hal_file = "./halLiftover_all_to_all/ref_based_small_chr21.hal"

    with Toil(options) as workflow:
        if not workflow.options.restart:
            #importing files:
            for asm, asm_file in assembly_files.items():
                assembly_files[asm] = workflow.importFile("file://" + os.path.abspath(asm_file))
            
            hal_file = workflow.importFile("file://" + os.path.abspath(options.hal_file))
                
            # output = workflow.start(Job.wrapJobFn(get_bases_unmapped, assembly_files, hal_file, options))
            output = workflow.start(Job.wrapJobFn(get_asm_mapping_depths, assembly_files, hal_file))
            
        else:
            output = workflow.restart()
        # write output
        print("output:", output)
        
        #todo: make output of get_asm_mapping_depths file-saveable.
        workflow.exportFile(output, 'file://' + os.path.abspath(options.output))

            


if __name__ == "__main__":
    main()