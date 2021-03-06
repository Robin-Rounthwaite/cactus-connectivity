#%%
"""
Built to exclude overlap in bedfiles with lots of overlap. (e.g. the halLiftover output)
"""
import collections as col
import operator

from types import SimpleNamespace
import ast

def get_mapping_coverage_points(alignment_bed):
    """
    Returns:
    all start and stop points of lines in the bedfile, sorted by contigs.
        key: contig_id, value: list[regions in tuple(point_value, start_bool) format].
        where start_bool is true if the point is a start of a region, and false if the point is a stop of the region.
    """
    # start-points and stop-points of each line in the bedfile. 
    # key: tuple(fasta_file, contig_id), value: list[regions in tuple(point_value, start_bool) format].
    mapping_coverage_points = col.defaultdict(list)

    # add all start and end points for regions that map well 
    with open(alignment_bed) as f:
        for line in f:
            # parse line in map_file:
            parsed = line.split("\t")
            
            contig_id = parsed[0]
            start = int(parsed[1])
            stop = int(parsed[2])

            # add these coordinates to mapping_coverage_points
            mapping_coverage_points[contig_id].append((start, True))
            mapping_coverage_points[contig_id].append((stop, False))
    return mapping_coverage_points

def get_mapping_coverage_coordinates(mapping_coverage_points):
    """
    Returns all the coords (defined by tuple(start,stop)) that are covered by at least one mapping in 
    mapping_coverage_points.
    """
    # mapping_coverage_coords is key: contig_id, value: list of coords: [(start, stop)]
    mapping_coverage_coords = col.defaultdict(list)
    for contig_id in mapping_coverage_points:
        contig_coverage_points = sorted(mapping_coverage_points[contig_id], key=operator.itemgetter(0, 1))
        # contig_coverage_points = sorted(mapping_coverage_points[contig_id], key=lambda point: point[0])
        open_points = 0
        current_region = [0, 0] # format (start, stop)
        for i in range(len(contig_coverage_points)):
            if open_points:
                # then we have at least one read overlapping this region.
                # expand the stop point of current_region
                current_region[1] = contig_coverage_points[i][0]
            if contig_coverage_points[i][1]:
                # if start_bool is true, the point represents a start of mapping
                open_points += 1
                if open_points == 1:
                    # that is, if we've found the starting point of a new current_region,
                    # so we should set the start of the current_region.
                    current_region[0] = contig_coverage_points[i][0]
            else:
                # if start_bool is not true, the point represents the end of a mapping.
                open_points -= 1
                if not open_points:
                    # if there's no more open_points in this region, then this is the 
                    # end of the current_region. Save current_region.
                    mapping_coverage_coords[contig_id].append(current_region.copy())
    return mapping_coverage_coords

def get_poor_mapping_coverage_coordinates(contig_lengths, mapping_coverage_coords, options):
    """
    mapping_coverage_coords is a dictionary of lists of coords in (start, stop) format.
    This function returns poor mapping coords, which is essentially the gaps between 
        those coords.
    example: mapping_coverage_coords{contig_1:[(3,5), (7, 9)]} would result in
                mapping_coverage_coords{contig_1:[(0,3), (5,7), (9, 11)]}, if contig_1 had a
                length of 11.
    variables:
        contig_lengths: A dictionary of the length of all the contigs in 
            {key: contig_id value: len(contig)} format.
        mapping_coverage_coords: a dictionary of lists of coords in 
            {key: contig_id, value:[(start, stop)]}
        sequence_context: an integer, representing the amount of sequence you would 
            want to expand each of the poor_mapping_coords by, to include context
            sequence for the poor mapping sequence. 
    """
    # poor_mapping_coords has key: contig_id, value list(tuple_of_positions(start, stop))
    poor_mapping_coords = col.defaultdict(list)
    for contig_id in contig_lengths:
        if contig_id in mapping_coverage_coords:
            if mapping_coverage_coords[contig_id][0][0] > 0:
                # if the first mapping region for the contig doesn't start at the start of
                # the contig, the first region is between the start of the contig and the 
                # start of the good_mapping_region.
                poor_mapping_stop = mapping_coverage_coords[contig_id][0][0] + options.sequence_context
                if poor_mapping_stop > contig_lengths[contig_id]:
                    poor_mapping_stop = contig_lengths[contig_id]
                if poor_mapping_stop - 0 >= options.minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((0, poor_mapping_stop))
                else:
                    pass
            for i in range(len(mapping_coverage_coords[contig_id]) - 1):
                # for every pair of mapping coords i and i + 1,
                # make a pair of (stop_from_ith_region, start_from_i+1th_region) to
                # represent the poor_mapping_coords. Include sequence_context as necessary.
                poor_mapping_start = mapping_coverage_coords[contig_id][i][1] - options.sequence_context
                if poor_mapping_start < 0:
                    poor_mapping_start = 0
                    
                poor_mapping_stop = mapping_coverage_coords[contig_id][i + 1][0] + options.sequence_context
                if poor_mapping_stop > contig_lengths[contig_id]:
                    poor_mapping_stop = contig_lengths[contig_id]

                if poor_mapping_stop - poor_mapping_start >= options.minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((poor_mapping_start, poor_mapping_stop))
                else:
                    pass
            if mapping_coverage_coords[contig_id][-1][1] < contig_lengths[contig_id]:
                # if the last mapping region for the contig stops before the end of
                # the contig, the last region is between the end of the mapping and the 
                # end of the contig.
                poor_mapping_start = mapping_coverage_coords[contig_id][-1][1] - options.sequence_context
                if poor_mapping_start < 0:
                    poor_mapping_start = 0
                if contig_lengths[contig_id] - poor_mapping_start >= options.minimum_size_remap: # implement size threshold.
                    poor_mapping_coords[contig_id].append((poor_mapping_start, contig_lengths[contig_id]))
                else:
                    pass

        else:
            # there isn't a good_mapping region for this contig. The full length of 
            # the contig belongs in poor_mapping_coords.
            poor_mapping_coords[contig_id].append((0, contig_lengths[contig_id]))
    return poor_mapping_coords

def get_seq_lengths_from_length_file(length_file):
    seq_lengths = dict()
    with open(length_file) as inf:
        for line in inf:
            parsed = line.split()
            seq_lengths[parsed[0]] = int(parsed[1])
    return seq_lengths

#%%
"""
General project idea:

Determining the quantity of "free" bases in the reference in the cactus graph. 
That is, how many bases in the reference aren't recorded as mapped to any other sequences
in the cactus graph?

Requires I merge together the mapping coverage points of liftovers from all three input
hg002 assemblies (where hg38 was the target, and the different assemblies the respective 
sources of each liftover.).
"""
def calc_free_bases(liftover_bed_files, contig_lengths, minimum_size_gap):
    mapping_coverage_points = col.defaultdict(list)
    for bedfile in liftover_bed_files:
        bed_mapping_coverage_points = get_mapping_coverage_points(bedfile)
        for key, values in bed_mapping_coverage_points.items():
            for value in values:
                mapping_coverage_points[key].append(value)
    mapping_coverage_coordinates = get_mapping_coverage_coordinates(mapping_coverage_points)

    options = SimpleNamespace()
    options.sequence_context = 0 #The nonzero default used in ref-based pipeline doesn't make sense in this context. 
    options.minimum_size_remap = minimum_size_gap

    poor_mapping_coverage_coordinates = get_poor_mapping_coverage_coordinates(contig_lengths, mapping_coverage_coordinates, options)

    seq_len_tot = sum(contig_lengths.values())
    print("Total length of reference:")
    print(seq_len_tot)
    print()

    unmapped_seq_len = int()
    for i in poor_mapping_coverage_coordinates.values():
        unmapped_seq_len += sum([int(j[1])-int(j[0]) for j in i])
    print("bases in target of liftover bed that are unaligned to the source sequence (with gaps between alignments <" + str(minimum_size_gap) + " bases in size still counted as aligned):")
    print(str(unmapped_seq_len) + " (" + str(round((unmapped_seq_len/seq_len_tot)*100, 2)) + "%)")


#%%
"""
For determining bases unmapped in hg38 reference. (hg38 is the target)
"""

file_prefix = "./hg002_data/ref_target_liftover_beds/"
# print(full_bed_files)

# print(bed_files)
bed_files = ["asm8a_liftover_hg38.bed", "asm9a_liftover_hg38.bed", "asm21_liftover_hg38.bed"]
# print(bed_files)
full_bed_files = [file_prefix + bed_file for bed_file in bed_files]
# print(bed_files)
# print(full_bed_files)

ref_contig_lengths = ast.literal_eval("{'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chrX': 156040895, 'chr8': 145138636, 'chr9': 138394717, 'chr11': 135086622, 'chr10': 133797422, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr20': 64444167, 'chr19': 58617616, 'chrY': 57227415, 'chr22': 50818468, 'chr21': 46709983, 'chr8_KZ208915v1_fix': 6367528, 'chr15_KI270905v1_alt': 5161414, 'chr15_KN538374v1_fix': 4998962, 'chr6_GL000256v2_alt': 4929269, 'chr6_GL000254v2_alt': 4827813, 'chr6_GL000251v2_alt': 4795265, 'chr6_GL000253v2_alt': 4677643, 'chr6_GL000250v2_alt': 4672374, 'chr6_GL000255v2_alt': 4606388, 'chr6_GL000252v2_alt': 4604811, 'chr17_KI270857v1_alt': 2877074, 'chr16_KI270853v1_alt': 2659700, 'chr15_KQ031389v1_alt': 2365364, 'chr16_KV880768v1_fix': 1927115, 'chr16_KI270728v1_random': 1872759, 'chr17_GL000258v2_alt': 1821992, 'chr5_GL339449v2_alt': 1612928, 'chr14_KI270847v1_alt': 1511111, 'chr17_KI270908v1_alt': 1423190, 'chr14_KI270846v1_alt': 1351393, 'chr5_KI270897v1_alt': 1144418, 'chr7_KI270803v1_alt': 1111570, 'chr19_GL949749v2_alt': 1091841, 'chr19_KI270938v1_alt': 1066800, 'chr19_GL949750v2_alt': 1066390, 'chr19_GL949748v2_alt': 1064304, 'chr12_KZ208916v1_fix': 1046838, 'chr19_GL949751v2_alt': 1002683, 'chr19_GL949746v1_alt': 987716, 'chr19_GL949752v1_alt': 987100, 'chr8_KI270821v1_alt': 985506, 'chr1_KI270763v1_alt': 911658, 'chr6_KI270801v1_alt': 870480, 'chr19_GL949753v2_alt': 796479, 'chr19_GL949747v2_alt': 729520, 'chr14_KZ208920v1_fix': 690932, 'chr7_KZ208913v1_alt': 680662, 'chr5_KV575244v1_fix': 673059, 'chr8_KI270822v1_alt': 624492, 'chr7_KZ208912v1_fix': 589656, 'chr4_GL000257v2_alt': 586476, 'chr12_KI270904v1_alt': 572349, 'chr4_KI270925v1_alt': 555799, 'chr1_KV880763v1_alt': 551020, 'chr12_KN538369v1_fix': 541038, 'chr2_KQ983256v1_alt': 535088, 'chr2_KQ031384v1_fix': 481245, 'chr16_KZ559113v1_fix': 480415, 'chr15_KI270852v1_alt': 478999, 'chr7_KV880765v1_fix': 468267, 'chr1_KQ031383v1_fix': 467143, 'chr1_KN538360v1_fix': 460100, 'chr3_KN196475v1_fix': 451168, 'chr15_KI270727v1_random': 448248, 'chr9_KI270823v1_alt': 439082, 'chr15_KI270850v1_alt': 430880, 'chr1_KI270759v1_alt': 425601, 'chr4_KV766193v1_alt': 420675, 'chr10_KN538367v1_fix': 420164, 'chr3_KN538364v1_fix': 415308, 'chr3_KV766192v1_fix': 411654, 'chr12_GL877876v1_alt': 408271, 'chr18_KQ090028v1_fix': 407387, 'chr19_KQ458386v1_fix': 405389, 'chrUn_KI270442v1': 392061, 'chr17_KI270862v1_alt': 391357, 'chr15_GL383555v2_alt': 388773, 'chr19_GL383573v1_alt': 385657, 'chr4_KI270896v1_alt': 378547, 'chr4_GL383528v1_alt': 376187, 'chr17_GL383563v3_alt': 375691, 'chr8_KI270810v1_alt': 374415, 'chr3_KQ031385v1_fix': 373699, 'chr19_KN196484v1_fix': 370917, 'chr1_GL383520v2_alt': 366580, 'chr2_KN538363v1_fix': 365499, 'chr5_KV575243v1_alt': 362221, 'chr13_KN538372v1_fix': 356766, 'chr1_KI270762v1_alt': 354444, 'chr1_KQ458383v1_alt': 349938, 'chr9_KN196479v1_fix': 330164, 'chr1_KZ208906v1_fix': 330031, 'chr15_KI270848v1_alt': 327382, 'chr17_KI270909v1_alt': 325800, 'chr14_KI270844v1_alt': 322166, 'chr6_KQ031387v1_fix': 320750, 'chr8_KI270900v1_alt': 318687, 'chr12_KQ759760v1_fix': 315610, 'chr10_GL383546v1_alt': 309802, 'chr13_KI270838v1_alt': 306913, 'chr3_KN196476v1_fix': 305979, 'chr8_KI270816v1_alt': 305841, 'chr1_KN538361v1_fix': 305542, 'chr11_KZ559108v1_fix': 305244, 'chr22_KI270879v1_alt': 304135, 'chr3_KZ559103v1_alt': 302885, 'chr11_KZ559110v1_alt': 301637, 'chr8_KI270813v1_alt': 300230, 'chr11_KI270831v1_alt': 296895, 'chr15_GL383554v1_alt': 296527, 'chr19_KV575249v1_alt': 293522, 'chr8_KI270811v1_alt': 292436, 'chr18_GL383567v1_alt': 289831, 'chrX_KI270880v1_alt': 284869, 'chr8_KI270812v1_alt': 282736, 'chr19_KI270921v1_alt': 282224, 'chr17_KV766196v1_fix': 281919, 'chr17_KI270729v1_random': 280839, 'chr11_KZ559109v1_fix': 279644, 'chr1_KQ983255v1_alt': 278659, 'chr17_JH159146v1_alt': 278131, 'chr10_KN196480v1_fix': 277797, 'chr17_KV766198v1_alt': 276292, 'chrX_KI270913v1_alt': 274009, 'chr6_KI270798v1_alt': 271782, 'chr7_KI270808v1_alt': 271455, 'chr6_KN196478v1_fix': 268330, 'chr16_KQ090027v1_alt': 267463, 'chr8_KV880767v1_fix': 265876, 'chr10_KQ090021v1_fix': 264545, 'chr22_KI270876v1_alt': 263666, 'chr15_KI270851v1_alt': 263054, 'chr22_KI270875v1_alt': 259914, 'chr1_KI270766v1_alt': 256271, 'chr19_KI270882v1_alt': 248807, 'chr3_KI270778v1_alt': 248252, 'chr17_KV766197v1_alt': 246895, 'chr6_KQ090016v1_fix': 245716, 'chr15_KI270849v1_alt': 244917, 'chr4_KI270786v1_alt': 244096, 'chr6_KZ208911v1_fix': 242796, 'chr19_KV575250v1_alt': 241058, 'chr12_KI270835v1_alt': 238139, 'chr4_KQ090015v1_alt': 236512, 'chr17_KI270858v1_alt': 235827, 'chr19_KI270867v1_alt': 233762, 'chr16_KI270855v1_alt': 232857, 'chr18_KZ559115v1_fix': 230843, 'chr4_KQ983257v1_fix': 230434, 'chr8_KI270926v1_alt': 229282, 'chr5_GL949742v1_alt': 226852, 'chr3_KI270780v1_alt': 224108, 'chr17_GL383565v1_alt': 223995, 'chr2_KI270774v1_alt': 223625, 'chr19_KV575256v1_alt': 223118, 'chr4_KI270790v1_alt': 220246, 'chr11_KI270927v1_alt': 218612, 'chr19_KI270932v1_alt': 215732, 'chr11_KI270903v1_alt': 214625, 'chr2_KI270894v1_alt': 214158, 'chr1_KQ458384v1_alt': 212205, 'chr12_KN196482v1_fix': 211377, 'chr14_GL000225v1_random': 211173, 'chrUn_KI270743v1': 210658, 'chr11_KI270832v1_alt': 210133, 'chr7_KI270805v1_alt': 209988, 'chrY_KZ208924v1_fix': 209722, 'chr4_GL000008v2_random': 209709, 'chr7_KI270809v1_alt': 209586, 'chr19_KI270887v1_alt': 209512, 'chr2_KN538362v1_fix': 208149, 'chr13_KN538371v1_fix': 206320, 'chr4_KI270789v1_alt': 205944, 'chr4_KQ983258v1_alt': 205407, 'chr3_KI270779v1_alt': 205312, 'chr19_KI270914v1_alt': 205194, 'chr18_KQ458385v1_alt': 205101, 'chr19_KI270886v1_alt': 204239, 'chr11_KI270829v1_alt': 204059, 'chr11_KN538368v1_alt': 203552, 'chr14_GL000009v2_random': 201709, 'chr21_GL383579v2_alt': 201197, 'chr11_JH159136v1_alt': 200998, 'chr19_KI270930v1_alt': 200773, 'chrUn_KI270747v1': 198735, 'chr18_GL383571v1_alt': 198278, 'chr19_KI270920v1_alt': 198005, 'chr3_KZ559102v1_alt': 197752, 'chr6_KI270797v1_alt': 197536, 'chr3_KI270935v1_alt': 197351, 'chr11_KQ759759v1_fix': 196940, 'chr17_KI270861v1_alt': 196688, 'chr15_KI270906v1_alt': 196384, 'chr5_KI270791v1_alt': 195710, 'chr3_KZ559105v1_alt': 195063, 'chr14_KI270722v1_random': 194050, 'chr16_GL383556v1_alt': 192462, 'chr13_KI270840v1_alt': 191684, 'chr14_GL000194v1_random': 191469, 'chr11_JH159137v1_alt': 191409, 'chr19_KI270917v1_alt': 190932, 'chr7_KI270899v1_alt': 190869, 'chr19_KI270923v1_alt': 189352, 'chr10_KI270825v1_alt': 188315, 'chr19_GL383576v1_alt': 188024, 'chrX_KV766199v1_alt': 188004, 'chr19_KI270922v1_alt': 187935, 'chrUn_KI270742v1': 186739, 'chr1_KN196472v1_fix': 186494, 'chr22_KI270878v1_alt': 186262, 'chr19_KI270929v1_alt': 186203, 'chr11_KI270826v1_alt': 186169, 'chr6_KB021644v2_alt': 185823, 'chr17_GL000205v2_random': 185591, 'chr10_KQ090020v1_alt': 185507, 'chr1_KI270765v1_alt': 185285, 'chr19_KI270916v1_alt': 184516, 'chr19_KI270890v1_alt': 184499, 'chr3_KI270784v1_alt': 184404, 'chr12_GL383551v1_alt': 184319, 'chr20_KI270870v1_alt': 183433, 'chrUn_GL000195v1': 182896, 'chr1_GL383518v1_alt': 182439, 'chr11_KQ090022v1_fix': 181958, 'chr22_KI270736v1_random': 181920, 'chr2_KZ208907v1_alt': 181658, 'chr10_KI270824v1_alt': 181496, 'chr11_KZ559111v1_alt': 181167, 'chr14_KI270845v1_alt': 180703, 'chr3_GL383526v1_alt': 180671, 'chr13_KI270839v1_alt': 180306, 'chr7_KQ031388v1_fix': 179932, 'chr22_KI270733v1_random': 179772, 'chrUn_GL000224v1': 179693, 'chr10_GL383545v1_alt': 179254, 'chrUn_GL000219v1': 179198, 'chr5_KI270792v1_alt': 179043, 'chr17_KI270860v1_alt': 178921, 'chr19_KV575252v1_alt': 178197, 'chr19_GL000209v2_alt': 177381, 'chr11_KI270830v1_alt': 177092, 'chr9_KI270719v1_random': 176845, 'chrUn_GL000216v2': 176608, 'chr22_KI270928v1_alt': 176103, 'chr1_KI270712v1_random': 176043, 'chr3_KZ208909v1_alt': 175849, 'chr6_KI270800v1_alt': 175808, 'chr1_KI270706v1_random': 175055, 'chr12_KZ208918v1_alt': 174808, 'chr22_KQ458388v1_alt': 174749, 'chr2_KI270776v1_alt': 174166, 'chr18_KI270912v1_alt': 174061, 'chr3_KI270777v1_alt': 173649, 'chr5_GL383531v1_alt': 173459, 'chr3_JH636055v2_alt': 173151, 'chr14_KI270725v1_random': 172810, 'chr5_KI270796v1_alt': 172708, 'chr7_KZ559106v1_alt': 172555, 'chr14_KZ208919v1_alt': 171798, 'chr9_GL383541v1_alt': 171286, 'chr19_KV575259v1_alt': 171263, 'chr19_KI270885v1_alt': 171027, 'chr19_KI270919v1_alt': 170701, 'chr19_KI270889v1_alt': 170698, 'chr19_KI270891v1_alt': 170680, 'chr19_KI270915v1_alt': 170665, 'chr19_KI270933v1_alt': 170537, 'chr19_KI270883v1_alt': 170399, 'chr19_GL383575v2_alt': 170222, 'chr19_KV575247v1_alt': 170206, 'chr19_KI270931v1_alt': 170148, 'chr12_GL383550v2_alt': 169178, 'chr16_KQ031390v1_alt': 169136, 'chr13_KI270841v1_alt': 169134, 'chrUn_KI270744v1': 168472, 'chr13_KQ090024v1_alt': 168146, 'chr19_KV575248v1_alt': 168131, 'chr18_KI270863v1_alt': 167999, 'chr18_GL383569v1_alt': 167950, 'chr12_GL877875v1_alt': 167313, 'chr21_KI270874v1_alt': 166743, 'chr19_KV575253v1_alt': 166713, 'chr3_KI270924v1_alt': 166540, 'chr1_KN196473v1_fix': 166200, 'chr1_KZ208904v1_alt': 166136, 'chr1_KI270761v1_alt': 165834, 'chr3_KQ031386v1_fix': 165718, 'chr3_KI270937v1_alt': 165607, 'chr8_KZ208914v1_fix': 165120, 'chr22_KI270734v1_random': 165050, 'chr18_GL383570v1_alt': 164789, 'chr5_KI270794v1_alt': 164558, 'chr4_GL383527v1_alt': 164536, 'chrUn_GL000213v1': 164239, 'chr3_KI270936v1_alt': 164170, 'chr3_KZ559101v1_alt': 164041, 'chr19_KV575246v1_alt': 163926, 'chr9_KQ090018v1_alt': 163882, 'chr4_KQ090014v1_alt': 163749, 'chr3_KI270934v1_alt': 163458, 'chr18_KZ559116v1_alt': 163186, 'chr9_GL383539v1_alt': 162988, 'chr3_KI270895v1_alt': 162896, 'chr22_GL383582v2_alt': 162811, 'chr3_KI270782v1_alt': 162429, 'chr1_KI270892v1_alt': 162212, 'chrUn_GL000220v1': 161802, 'chr2_KI270767v1_alt': 161578, 'chr2_KI270715v1_random': 161471, 'chr2_KI270893v1_alt': 161218, 'chrUn_GL000218v1': 161147, 'chr19_KV575255v1_alt': 161095, 'chr18_GL383572v1_alt': 159547, 'chr19_KV575251v1_alt': 159285, 'chr8_KI270817v1_alt': 158983, 'chr4_KI270788v1_alt': 158965, 'chrUn_KI270749v1': 158759, 'chr7_KI270806v1_alt': 158166, 'chr7_KI270804v1_alt': 157952, 'chr18_KI270911v1_alt': 157710, 'chrUn_KI270741v1': 157432, 'chr17_KI270910v1_alt': 157099, 'chr19_KI270884v1_alt': 157053, 'chr8_KV880766v1_fix': 156998, 'chr19_KV575258v1_alt': 156965, 'chr22_KN196485v1_alt': 156562, 'chr22_KQ458387v1_alt': 155930, 'chr19_GL383574v1_alt': 155864, 'chr19_KI270888v1_alt': 155532, 'chr3_GL000221v1_random': 155397, 'chr17_KV575245v1_fix': 154723, 'chr11_GL383547v1_alt': 154407, 'chr12_KZ559112v1_alt': 154139, 'chr2_KI270716v1_random': 153799, 'chr22_KN196486v1_alt': 153027, 'chr12_GL383553v2_alt': 152874, 'chr6_KI270799v1_alt': 152148, 'chr22_KI270731v1_random': 150754, 'chrUn_KI270751v1': 150742, 'chrUn_KI270750v1': 148850, 'chr13_KN538373v1_fix': 148762, 'chr19_KV575260v1_alt': 145691, 'chr8_KI270818v1_alt': 145606, 'chr22_KQ759761v1_alt': 145162, 'chrX_KI270881v1_alt': 144206, 'chr21_KI270873v1_alt': 143900, 'chr2_GL383521v1_alt': 143390, 'chr7_KV880764v1_fix': 142129, 'chr8_KI270814v1_alt': 141812, 'chr1_KQ458382v1_alt': 141019, 'chr11_KV766195v1_fix': 140877, 'chr2_KZ208908v1_alt': 140361, 'chr1_KZ208905v1_alt': 140355, 'chr6_KV766194v1_fix': 139427, 'chr5_KN196477v1_alt': 139087, 'chr12_GL383552v1_alt': 138655, 'chrUn_KI270519v1': 138126, 'chr2_KI270775v1_alt': 138019, 'chr17_KI270907v1_alt': 137721, 'chrUn_GL000214v1': 137718, 'chr8_KI270901v1_alt': 136959, 'chr2_KI270770v1_alt': 136240, 'chr5_KZ208910v1_alt': 135987, 'chr16_KI270854v1_alt': 134193, 'chr9_KQ090019v1_alt': 134099, 'chr8_KI270819v1_alt': 133535, 'chr17_GL383564v2_alt': 133151, 'chr2_KI270772v1_alt': 133041, 'chr8_KI270815v1_alt': 132244, 'chr5_KI270795v1_alt': 131892, 'chr5_KI270898v1_alt': 130957, 'chr20_GL383577v2_alt': 128386, 'chr1_KI270708v1_random': 127682, 'chr7_KI270807v1_alt': 126434, 'chr5_KI270793v1_alt': 126136, 'chr6_GL383533v1_alt': 124736, 'chr2_GL383522v1_alt': 123821, 'chr13_KQ090025v1_alt': 123480, 'chr19_KI270918v1_alt': 123111, 'chr1_KN196474v1_fix': 122022, 'chr12_GL383549v1_alt': 120804, 'chr2_KI270769v1_alt': 120616, 'chr4_KI270785v1_alt': 119912, 'chr12_KI270834v1_alt': 119498, 'chr7_GL383534v2_alt': 119183, 'chr20_KI270869v1_alt': 118774, 'chr17_KZ559114v1_alt': 116753, 'chr21_GL383581v2_alt': 116689, 'chr3_KI270781v1_alt': 113034, 'chr17_KI270730v1_random': 112551, 'chrUn_KI270438v1': 112505, 'chr4_KI270787v1_alt': 111943, 'chr18_KI270864v1_alt': 111737, 'chr2_KI270771v1_alt': 110395, 'chr1_GL383519v1_alt': 110268, 'chr2_KI270768v1_alt': 110099, 'chr1_KI270760v1_alt': 109528, 'chr12_KQ090023v1_alt': 109323, 'chr3_KI270783v1_alt': 109187, 'chr11_KN196481v1_fix': 108875, 'chr17_KI270859v1_alt': 108763, 'chr11_KI270902v1_alt': 106711, 'chr3_KZ559104v1_fix': 105527, 'chr18_GL383568v1_alt': 104552, 'chr22_KI270737v1_random': 103838, 'chr13_KI270843v1_alt': 103832, 'chr8_KZ559107v1_alt': 103072, 'chr22_KI270877v1_alt': 101331, 'chr5_GL383530v1_alt': 101241, 'chrY_KN196487v1_fix': 101150, 'chr22_KQ759762v1_fix': 101037, 'chr19_KV575257v1_alt': 100553, 'chr11_KI270721v1_random': 100316, 'chr19_KV575254v1_alt': 99845, 'chr22_KI270738v1_random': 99375, 'chr22_GL383583v2_alt': 96924, 'chr2_GL582966v2_alt': 96131, 'chrUn_KI270748v1': 93321, 'chr18_KZ208922v1_fix': 93070, 'chrUn_KI270435v1': 92983, 'chr5_GL000208v1_random': 92689, 'chrUn_KI270538v1': 91309, 'chr4_KQ090013v1_alt': 90922, 'chr17_GL383566v1_alt': 90219, 'chr16_GL383557v1_alt': 89672, 'chr17_JH159148v1_alt': 88070, 'chr12_KN538370v1_fix': 86533, 'chr10_KN538366v1_fix': 85284, 'chr5_GL383532v1_alt': 82728, 'chr21_KI270872v1_alt': 82692, 'chr6_KQ090017v1_alt': 82315, 'chrUn_KI270756v1': 79590, 'chr16_KZ208921v1_alt': 78609, 'chr6_KI270758v1_alt': 76752, 'chr12_KI270833v1_alt': 76061, 'chr6_KI270802v1_alt': 75005, 'chr21_GL383580v2_alt': 74653, 'chr22_KB663609v1_alt': 74013, 'chr22_KI270739v1_random': 73985, 'chr9_GL383540v1_alt': 71551, 'chrUn_KI270757v1': 71251, 'chr2_KI270773v1_alt': 70887, 'chr17_JH159147v1_alt': 70345, 'chr11_KI270827v1_alt': 67707, 'chr1_KI270709v1_random': 66860, 'chrUn_KI270746v1': 66486, 'chr12_KZ208917v1_fix': 64689, 'chr16_KI270856v1_alt': 63982, 'chr21_GL383578v2_alt': 63917, 'chrUn_KI270753v1': 62944, 'chr19_KI270868v1_alt': 61734, 'chr9_GL383542v1_alt': 60032, 'chr16_KQ090026v1_alt': 59016, 'chr20_KI270871v1_alt': 58661, 'chr12_KI270836v1_alt': 56134, 'chr19_KI270865v1_alt': 52969, 'chr1_KI270764v1_alt': 50258, 'chrY_KZ208923v1_fix': 48370, 'chr1_KZ559100v1_fix': 44955, 'chrUn_KI270589v1': 44474, 'chr14_KI270726v1_random': 43739, 'chr19_KI270866v1_alt': 43156, 'chr22_KI270735v1_random': 42811, 'chr1_KI270711v1_random': 42210, 'chrUn_KI270745v1': 41891, 'chr1_KI270714v1_random': 41717, 'chr22_KI270732v1_random': 41543, 'chr1_KI270713v1_random': 40745, 'chrUn_KI270754v1': 40191, 'chr1_KI270710v1_random': 40176, 'chr12_KI270837v1_alt': 40090, 'chr9_KI270717v1_random': 40062, 'chr14_KI270724v1_random': 39555, 'chr9_KI270720v1_random': 39050, 'chr14_KI270723v1_random': 38115, 'chr9_KI270718v1_random': 38054, 'chrUn_KI270317v1': 37690, 'chr13_KI270842v1_alt': 37287, 'chrY_KI270740v1_random': 37240, 'chrUn_KI270755v1': 36723, 'chr8_KI270820v1_alt': 36640, 'chr13_KN196483v1_fix': 35455, 'chr1_KI270707v1_random': 32032, 'chrUn_KI270579v1': 31033, 'chrUn_KI270752v1': 27745, 'chrUn_KI270512v1': 22689, 'chrUn_KI270322v1': 21476, 'chrM': 16569, 'chrUn_GL000226v1': 15008, 'chr10_KN538365v1_fix': 14347, 'chrUn_KI270311v1': 12399, 'chrUn_KI270366v1': 8320, 'chrUn_KI270511v1': 8127, 'chrUn_KI270448v1': 7992, 'chrUn_KI270521v1': 7642, 'chrUn_KI270581v1': 7046, 'chrUn_KI270582v1': 6504, 'chrUn_KI270515v1': 6361, 'chrUn_KI270588v1': 6158, 'chrUn_KI270591v1': 5796, 'chrUn_KI270522v1': 5674, 'chrUn_KI270507v1': 5353, 'chrUn_KI270590v1': 4685, 'chrUn_KI270584v1': 4513, 'chrUn_KI270320v1': 4416, 'chrUn_KI270382v1': 4215, 'chrUn_KI270468v1': 4055, 'chrUn_KI270467v1': 3920, 'chrUn_KI270362v1': 3530, 'chrUn_KI270517v1': 3253, 'chrUn_KI270593v1': 3041, 'chrUn_KI270528v1': 2983, 'chrUn_KI270587v1': 2969, 'chrUn_KI270364v1': 2855, 'chrUn_KI270371v1': 2805, 'chrUn_KI270333v1': 2699, 'chrUn_KI270374v1': 2656, 'chrUn_KI270411v1': 2646, 'chrUn_KI270414v1': 2489, 'chrUn_KI270510v1': 2415, 'chrUn_KI270390v1': 2387, 'chrUn_KI270375v1': 2378, 'chrUn_KI270420v1': 2321, 'chrUn_KI270509v1': 2318, 'chrUn_KI270315v1': 2276, 'chrUn_KI270302v1': 2274, 'chrUn_KI270518v1': 2186, 'chrUn_KI270530v1': 2168, 'chrUn_KI270304v1': 2165, 'chrUn_KI270418v1': 2145, 'chrUn_KI270424v1': 2140, 'chrUn_KI270417v1': 2043, 'chrUn_KI270508v1': 1951, 'chrUn_KI270303v1': 1942, 'chrUn_KI270381v1': 1930, 'chrUn_KI270529v1': 1899, 'chrUn_KI270425v1': 1884, 'chrUn_KI270396v1': 1880, 'chrUn_KI270363v1': 1803, 'chrUn_KI270386v1': 1788, 'chrUn_KI270465v1': 1774, 'chrUn_KI270383v1': 1750, 'chrUn_KI270384v1': 1658, 'chrUn_KI270330v1': 1652, 'chrUn_KI270372v1': 1650, 'chrUn_KI270548v1': 1599, 'chrUn_KI270580v1': 1553, 'chrUn_KI270387v1': 1537, 'chrUn_KI270391v1': 1484, 'chrUn_KI270305v1': 1472, 'chrUn_KI270373v1': 1451, 'chrUn_KI270422v1': 1445, 'chrUn_KI270316v1': 1444, 'chrUn_KI270338v1': 1428, 'chrUn_KI270340v1': 1428, 'chrUn_KI270583v1': 1400, 'chrUn_KI270334v1': 1368, 'chrUn_KI270429v1': 1361, 'chrUn_KI270393v1': 1308, 'chrUn_KI270516v1': 1300, 'chrUn_KI270389v1': 1298, 'chrUn_KI270466v1': 1233, 'chrUn_KI270388v1': 1216, 'chrUn_KI270544v1': 1202, 'chrUn_KI270310v1': 1201, 'chrUn_KI270412v1': 1179, 'chrUn_KI270395v1': 1143, 'chrUn_KI270376v1': 1136, 'chrUn_KI270337v1': 1121, 'chrUn_KI270335v1': 1048, 'chrUn_KI270378v1': 1048, 'chrUn_KI270379v1': 1045, 'chrUn_KI270329v1': 1040, 'chrUn_KI270419v1': 1029, 'chrUn_KI270336v1': 1026, 'chrUn_KI270312v1': 998, 'chrUn_KI270539v1': 993, 'chrUn_KI270385v1': 990, 'chrUn_KI270423v1': 981, 'chrUn_KI270392v1': 971, 'chrUn_KI270394v1': 970}")
minimum_size_remap = 0
# minimum_size_remap = 100
# minimum_size_remap = 10000

calc_free_bases(full_bed_files, ref_contig_lengths, minimum_size_remap)

"""Output:
Total length of reference (calculated as seq_len_tot):
3257347282


Length of sequence "free" (parts of reference with no mapping) with...
<10k base gaps filled/skipped:
296438393 (9.1%)
<100 base gaps filled/skipped:
346356570 (10.6%)
no base gaps filled (count all "free" regions of the reference, even if 1 bp in length):
348672674 (10.7%)
"""

#%%
"""
For determining bases unmapped to reference in each of three assemblies.
"""
bed_file_prefix = "./hg002_data/asm_target_liftover_beds/"
length_file_prefix = "./hg002_data/asm_lengths/"
all_bed_files = {"asm8a": "hg38_liftover_asm8a.bed", "asm9a": "hg38_liftover_asm9a.bed", "asm21": "hg38_liftover_asm21.bed"}
all_asm_lengths = {"asm8a": "asm8a.lengths", "asm9a": "asm9a.lengths", "asm21": "asm21.lengths"}

for asm in all_bed_files:
    # only need one bed_file because we're just looking at how each assembly(target) is 
    # aligned to the reference.
    bed_file = all_bed_files[asm]
    bed_files = [bed_file_prefix + bed_file]

    length_file = length_file_prefix + all_asm_lengths[asm]

    minimum_size_gap = 0
    assembly_seq_lengths = get_seq_lengths_from_length_file(length_file)
    
    print("Calculating bases unmapped in cactus graph for " + asm + ":")
    calc_free_bases(bed_files, assembly_seq_lengths, minimum_size_gap)
    print()
    print()

    

#%%
"""
For determining bases unmapped to all other sequences (reference + assemblies) in the three assemblies. 

(requires doing further liftovers between each pair of sequences. Do in k8s job?)
"""



#%%





























#%%
# len_tot = sum(contig_lengths.values())
# # len_tot_manual = sum(contig_lengths_manual_calc.values())
# # print(len_tot, len_tot_manual)
# print(len_tot)

# #%%
# for key in poor_mapping_coverage_coordinates:
#     prev_end = int()
#     for pair in poor_mapping_coverage_coordinates[key]:
#         if pair[0] < prev_end:
#             print("problem in ", key)
#         if pair[1] > contig_lengths[key]:
#             print("points past contig lengths")
#         prev_end = pair[1]
        
# print(type(mapping_coverage_coordinates))

# #%%

# bedfile = "asm9a_liftover.bed"
# mapping_coverage_points = get_mapping_coverage_points(bedfile)
# asm9a_liftover_dipcalls_merged = get_mapping_coverage_coordinates(mapping_coverage_points)
# print(asm9a_liftover_dipcalls_merged)
# #%%
# # print(len(asm9a_liftover_dipcalls_merged.keys()))

# interval_sum = int()
# for contig, coords in asm9a_liftover_dipcalls_merged.items():
#     for coord in coords:
#         interval_sum += coord[1] - coord[0]

# #%%
# print(interval_sum)
# #%%
# # bedfile = "asm9ab.dip.bed"
# bedfile = "asm9a_liftover.bed"


# interval_sum = int()
# prev_interval_end = int()
# prev_interval_contig = str()
# for line in open(bedfile):
#     # print("hi")
#     parsed = line.split()
#     # if prev_interval_end > int(parsed[1]) and prev_interval_contig == parsed[0]:
#     #     print("----------------we have overlap! Prev interval ended at:", prev_interval_end, "current interval is ", parsed[0], parsed[1], parsed[2])

#     prev_interval_end = int(parsed[2])
#     prev_interval_contig = parsed[0]
#     interval_sum += int(parsed[2]) - int(parsed[1])

#     #Naive accounting for overlap: (doesnt work)
#     # if prev_interval_end > int(parsed[1]) and prev_interval_contig == parsed[0]:
#     #     if int(parsed[2]) > prev_interval_end:
#     #         interval_sum += int(parsed[2]) - prev_interval_end
#     #     else:
#     #         print("skipping interval altogether:", prev_interval_contig, prev_interval_end, parsed[0], parsed[1], parsed[2])
#     # else:
#     #     interval_sum += int(parsed[2]) - int(parsed[1]) 
        

# # %%
# print(interval_sum)


# """
# output:
# print(interval_sum) for asm9ab.dip.bed:
# 2815317002
# These are the regions of the 

# for asm9a_liftover.txt, without accounting for the overlap (which makes many of the bases far multiple-counted)
# 2800677903
# tried to naively account for overlap, got
# 0
# Used get_mapping_coverage_points trick, got:
# 2714090357
# """

# # %%
# options.sequence_context = 3
# print(options.sequence_context)

# #%%

# sn = SimpleNamespace()
# sn.a = test
# sn.a
# #%%

# %%
print("hi")

# %%
