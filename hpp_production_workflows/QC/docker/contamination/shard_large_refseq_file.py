#!/usr/bin/env python3

import sys
from Bio import SeqIO

def main():

	# required arguments
	if len(sys.argv) < 3:
		exit('Wrong number of arguments')
	fasta_input_file = sys.argv[1]
	fasta_output_prefix = sys.argv[2]

	# optional output size param
	desired_output_size = 2**30 if len(sys.argv) < 4 else int(sys.argv[3])

	# prep for IO
	fasta_output_handle = None
	fasta_input_handle = None
	current_output_size = 0
	current_output_index = 0

	try:
		# open sequence file and iterate over entries
		fasta_input_handle = open(fasta_input_file, 'r')
		for record in SeqIO.parse(fasta_input_handle, "fasta"):

			# maybe create a new file
			if fasta_output_handle is None or current_output_size > desired_output_size:
				if fasta_output_handle is not None:
					fasta_output_handle.close()
				filename = "{}_{}.fa".format(fasta_output_prefix, current_output_index)
				print("Writing to file {}".format(filename))
				fasta_output_handle = open(filename, 'w')
				current_output_index += 1
				current_output_size = 0

			# track seq length
			length = len(record.seq)
			current_output_size += length

			# write
			SeqIO.write([record], fasta_output_handle, 'fasta')


	except:
		if fasta_output_handle is not None:
			fasta_output_handle.close()
		if fasta_input_handle is not None:
			fasta_input_handle.close()

	print("Fin.")

	
if __name__ == "__main__":
	main()

