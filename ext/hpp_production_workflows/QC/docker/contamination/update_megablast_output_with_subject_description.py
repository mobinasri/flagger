#!/usr/bin/env python3

import sys

def main():
	if len(sys.argv) not in (4, 5):
		exit('Wrong number of arguments')

	reference_fasta = sys.argv[1]
	megablast_input_file = sys.argv[2]
	megablast_output_file = sys.argv[3]
	ancillary_info = None if len(sys.argv) == 4 else sys.argv[4]

	# build map
	id_map = dict()
	with open(reference_fasta, 'r') as reference_fasta_input:
		for line in reference_fasta_input:
			if line.startswith(">"):
				sequence_name = line[1:].split(" ")[0]
				sequence_metadata = " ".join(line.strip()[1:].split(" ")[1:])
				id_map[sequence_name] = sequence_metadata

	header_count = 0
	empty_count = 0
	found_subject_count = 0
	missing_subject_count = 0
	malformed_count = 0
	with open(megablast_input_file, 'r') as fin, open(megablast_output_file, 'w') as fout:
		for line in fin:
			line = line.strip()
			if len(line) == 0:
				empty_count += 1
			elif line.startswith("#"):
				line = line + "\tsubject_metadata"
				header_count += 1
			else:
				parts = line.split('\t')
				if len(parts) < 2:
					malformed_count += 1
				elif parts[1] in id_map:
					if ancillary_info is not None:
						line = line + "\t({}) {}".format(ancillary_info, id_map[parts[1]])
					else:
						line = line + "\t{}".format(id_map[parts[1]])
					found_subject_count += 1
				else:
					if ancillary_info is not None:
						line = line + "\t({}) N/A".format(ancillary_info)
					else:
						line = line + "\tN/A"
					missing_subject_count += 1
			fout.write(line + "\n")

	print("Updated megablast output with reference sequence description:", file=sys.stderr)
	print("\tHeader:    {}".format(header_count), file=sys.stderr)
	print("\tFound:     {}".format(found_subject_count), file=sys.stderr)
	print("\tMissing:   {}".format(missing_subject_count), file=sys.stderr)
	print("\tEmpty:     {}".format(empty_count), file=sys.stderr)
	print("\tMalformed: {}".format(malformed_count), file=sys.stderr)


if __name__ == "__main__":
	main()

