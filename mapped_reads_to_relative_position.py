#!/usr/bin/env python2.7

# M.D. Warman
# May 15, 2017

# This script takes a list of stop codon positions and a .bed file and outputs
# the cumulative read count (based on read start position) for each nt relative
# to the stop codon. Output is a two column tab-delimited file (binned_reads.txt):
#
# <1> position relative to stop codon
# <2> number of reads in this bin
#
# Usage:
# mapped_reads_to_relative_position.py [stop_codon_positions.txt] [reads.bed]

import sys
import io

if len(sys.argv) == 1:
    print("\nNo files supplied.\n\nUsage: "
     "mapped_reads_to_relative_position.py [stop_codon_positions.txt] "
     "[reads.bed]\n")

# These variables assign the distance upstream and downstream from the stop
# codon where reads will be considered part of that gene.
window_upstream_distance = 1000
window_downstream_distance = 1000

# Setting up some lists to hold the information associated with the 
# stop_codon_positions.txt file:
stops_file_chr = list()
stops_file_stop = list()
stops_file_strand = list()
window_start = list()
window_end = list()

# This list will contain counts for reads that map to every position within the
# region of interest relative to the stop codons. I'll just start by filling
# all positions with 0, not sure if that's necessary.
master_count = list()
master_position = list()
for x in range(2000):
    master_count.append(0)
for x in range(-1000, 1000):
    master_position.append(x)

# Reading in the stop codon positions file and setting up the window.
fandle_stops = io.open(sys.argv[1], "rU")
for line in fandle_stops:
    linestripped = line.strip()
    line_list = linestripped.split("\t")
    stops_file_chr.append(line_list[0])
    stops_file_stop.append(line_list[2])
    stops_file_strand.append(line_list[3])
    window_start.append(int(line_list[2]) - window_upstream_distance)
    window_end.append(int(line_list[2]) + window_upstream_distance)
fandle_stops.close()

# Here the reads file will be opened. For each read, the script will move
# through the stop codon position lists searching for a window that contains
# it. If a proper window is found, the distance from that stop codon will be
# calculated, and the read will be added to that bin in the master list.
fandle_reads_in = io.open(sys.argv[2], "rU")
line_counter = 0
old_read_pos = -1
for line in fandle_reads_in:
    linestripped = line.strip()
    line_list = linestripped.split("\t")
    read_chr = line_list[0]
    read_pos = int(line_list[1])
    # Doesn't do the whole search over again if the read position is the same 
    # previous read position
    if read_pos == old_read_pos:
        current_count = master_count[position_in_master_count]
        new_count = current_count + 1
        master_count[position_in_master_count] = new_count
    elif  line_list[5] == "+":         # Plus stranded reads
        # Searches for the proper window
        for pos in range(len(stops_file_stop)):
            current_chr = stops_file_chr[pos]
            current_stop = int(stops_file_stop[pos])
            current_strand = stops_file_strand[pos]
            current_win_start = int(window_start[pos])
            current_win_end = int(window_end[pos])
            if (read_chr == current_chr) and (read_pos < current_win_end) \
             and (read_pos > current_win_start):
                # Finds relative position and updates that position's count by
                # adding one. In the future it will add the normalized value
                # for that gene.
                old_read_pos = read_pos
                relative_position = read_pos - current_stop 
                position_in_master_count = relative_position + \
                 window_upstream_distance
                current_count = master_count[position_in_master_count]
                new_count = current_count + 1
                master_count[position_in_master_count] = new_count
                # Allows the script to stop searching once an appropriate window is
                # found.
                break
            #else:
            #   print("No match")
            #   print("Current chr: " + str(current_chr))
            #   print("Read chr: " + str(read_chr))
            #   print("Current position: " + str(read_pos))
            #   print("Current window: " + str(current_win_start) + \
            #    "\t" + str(current_win_end) + "\n\n")
fandle_reads_in.close()

# Writing the output.
out_fandle = io.open("binned_reads.txt", "wb")
for x in range(len(master_position)):
    out_fandle.write(str(master_position[x]) + "\t")
    out_fandle.write(str(master_count[x]) + "\n")
out_fandle.close()









