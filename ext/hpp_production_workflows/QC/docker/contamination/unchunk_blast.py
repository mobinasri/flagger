#!/usr/bin/env python3

import sys

CONTIG_IDX=0
START_IDX=6
END_IDX=7

def main():

	if len(sys.argv) < 3:
		exit('Wrong number of arguments')

	blast_input_file = sys.argv[1]
	blast_output_file = sys.argv[2]

	threshold_length = 1000000 if len(sys.argv) < 4 else int(sys.argv[3])

	with open(blast_input_file, 'r') as fin, open(blast_output_file, 'w') as fout:
		for i, line in enumerate(fin):
			if line.startswith("#"):
				fout.write(line)
				continue

			parts = line.split("\t")

			# hits have 4 parts, we don't care about misses
			if len(parts) < 4: continue
			chunked_contig = parts[CONTIG_IDX]
			chunked_start = int(parts[START_IDX])
			chunked_end = int(parts[END_IDX])

			# get chunking parts
			chunked_contig_parts = chunked_contig.split(".")
			if len(parts) <= 1 or not chunked_contig_parts[-1].startswith("chunk_"):
				raise Exception("Error parsing chunked annotation at line {} in {}: {}".format(i, blast_input_file, line))
			chunk_idx = int(chunked_contig_parts[-1].lstrip("chunk_"))

			# write it out
			parts[CONTIG_IDX] = ".".join(chunked_contig_parts[0:-1])
			parts[START_IDX] = str(chunked_start + threshold_length * chunk_idx)
			parts[END_IDX] = str(chunked_end + threshold_length * chunk_idx)
			fout.write("\t".join(parts))


if __name__ == "__main__":
	main()

