#!/usr/bin/env python3

import sys
from Bio import SeqIO

def main():

	if len(sys.argv) < 3:
		exit('Wrong number of arguments')

	fasta_input_file = sys.argv[1]
	fasta_output_file = sys.argv[2]

	threshold_length = 1000000 if len(sys.argv) < 4 else int(sys.argv[3])
	overlap_length = 10000 if len(sys.argv) < 5 else int(sys.argv[4])

	minimum_record_size = 11


	fasta_output_handle = open(fasta_output_file, 'w')

	with open(fasta_input_file, 'r') as fasta_input_handle:
		for record in SeqIO.parse(fasta_input_handle, "fasta"):

			if len(record) >= minimum_record_size:
				records_to_write = []

				slice_count = 0
				while (slice_count * threshold_length) < len(record) - (threshold_length+overlap_length):
					record_slice = record[(slice_count*threshold_length):((slice_count+1)*threshold_length + overlap_length)]
					record_slice.id += '.chunk_{}'.format(str(slice_count+1))

					record_slice.description = ''
					records_to_write.append(record_slice)
					slice_count += 1

				final_record_slice = record[(slice_count*threshold_length):]
				final_record_slice.id += '.chunk_' + str(slice_count+1)
				final_record_slice.description = ''

				records_to_write.append(final_record_slice)

				SeqIO.write(records_to_write, fasta_output_handle, 'fasta')

	fasta_output_handle.close()

if __name__ == "__main__":
	main()

