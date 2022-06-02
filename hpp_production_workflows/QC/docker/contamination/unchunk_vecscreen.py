#!/usr/bin/env python3

import sys

def main():

	if len(sys.argv) < 3:
		exit('Wrong number of arguments')

	vecscreen_input_file = sys.argv[1]
	vecscreen_output_file = sys.argv[2]

	threshold_length = 1000000 if len(sys.argv) < 4 else int(sys.argv[3])

	with open(vecscreen_input_file, 'r') as fin, open(vecscreen_output_file, 'w') as fout:
		for i, line in enumerate(fin):
			parts = line.split()

			# hits have 4 parts, we don't care about misses
			if len(parts) < 4: continue
			vecscreen_annotation = parts[0]
			chunked_contig = parts[1]
			chunked_start = int(parts[2])
			chunked_end = int(parts[3])

			# get chunking parts
			chunked_contig_parts = chunked_contig.split(".")
			if len(parts) <= 1 or not chunked_contig_parts[-1].startswith("chunk_"):
				raise Exception("Error parsing chunked annotation at line {} in {}: {}".format(i,vecscreen_input_file, line))
			chunk_idx = int(chunked_contig_parts[-1].lstrip("chunk_"))
			contig = ".".join(chunked_contig_parts[0:-1])

			# write it out
			start = chunked_start + threshold_length * chunk_idx
			end = chunked_end + threshold_length * chunk_idx
			fout.write("{}\t{}\t{}\t{}\n".format(contig, start, end, vecscreen_annotation))


if __name__ == "__main__":
	main()

