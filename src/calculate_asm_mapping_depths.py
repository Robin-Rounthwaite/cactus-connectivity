from src import calculate_bases_unmapped

import collections as col
import operator

def get_mapping_depths(job, mapping_coverage_points, contig_lengths):
    """
    Based on get_mapping_coverage_coordinates algorithm.
    Returns the number of bases covered at each depth level.
    """
    # mapping_depths is key: depth_level (int); value:bases_covered_at_depth_level
    # it measure the number of bases involved in an alignment, segregated by the number of
    # times each base is involved in an alignment.
    mapping_depths = col.defaultdict(int)
    
    for contig_id in mapping_coverage_points:
        liftover_coverage_points = sorted(mapping_coverage_points[contig_id], key=operator.itemgetter(0, 1))

        depth_coverage = 0
        last_base = 0 # The "previous" bases we measured ended at the beginning of the sequence. So, last_base is at 0.
        debug_1_if = int()
        debug_2_if = int()
        for point in liftover_coverage_points:
            # first, check to see if we've actually moved anywhere along the seq:
            if point[0] == last_base:
                if point[1]:
                    depth_coverage += 1
                else:
                    depth_coverage -= 1
                debug_1_if += 1
            else:
                # We've reached a point that's further along the seq.
                # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++mapping_depths", mapping_depths, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
                # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++depth_coverage", depth_coverage, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
                # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++last_base", last_base, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
                # print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++point", point, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
                mapping_depths[depth_coverage] += (point[0] - last_base)
                if point[1]:
                    depth_coverage += 1
                else:
                    depth_coverage -= 1
                last_base = point[0]
                debug_2_if += 1


        # check to make sure that any sequence after the last coverage point is included
        # in the mapping_depths dict at depth = 0.
        # if last_base + 1 < contig_lengths[contig_id]:
        if last_base < contig_lengths[contig_id]:
            # mapping_depths[0] += contig_lengths[contig_id] - last_base + 1
            mapping_depths[0] += contig_lengths[contig_id] - last_base

    return (mapping_depths, debug_1_if, debug_2_if)

def calculate_mapping_depths(job, liftover_bed_files, contig_lengths):
    # perform a separate calculation of intervals unmapped in each liftover_bed.
    # Then, add all the intervals into a single list, sorted by first digit, and then 
    # second digit.
    leader = job.addChildJobFn(calculate_bases_unmapped.empty)
    
    mapping_coverage_points = list()
    for bedfile in liftover_bed_files:
        mapping_coverage_points.append(leader.addChildJobFn(calculate_bases_unmapped.get_mapping_coverage_points, bedfile).rv())
    coverage_points_jobs = leader.encapsulate()

    merged_mapping_coverage_points = coverage_points_jobs.addChildJobFn(calculate_bases_unmapped.merge_mapping_coverage_points, mapping_coverage_points).rv()
    merging_jobs = coverage_points_jobs.encapsulate()

    mapping_depths = merging_jobs.addChildJobFn(get_mapping_depths, merged_mapping_coverage_points, contig_lengths).rv()
    mapping_depths_job = merging_jobs.encapsulate()

    return mapping_depths

def calculate_all_mapping_depths(job, liftovers, contig_lengths):
    #todo: implement minimum_size_gap, similar to in calculate_bases_unmapped?
    # mapping_depths has key: assembly_id value:list(bases_unmapped, bases_mapped_once, bases_mapped_twice... etc.)
    mapping_depths = dict()
    for target_assembly, source_assembly_liftovers in liftovers.items():
        mapping_depths[target_assembly] = job.addChildJobFn(calculate_mapping_depths, list(source_assembly_liftovers.values()), contig_lengths[target_assembly]).rv()
    return mapping_depths


"""    
Scratch:
# Then, walk through the list of intervals, saving each base as a separate value in a dict,
# until we've walked past all intervals that include that given base. Once we've counted the amount of coverage
# on that base, we pop the base and 


    # NOTE: if a base is counted as covered at 3x depth coverage, then it will also be counted as
    # 2x and 1x depth coverage as well. If undesirable, this is easily resolved with subtraction.

"""
