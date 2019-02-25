import gzip, argparse

def main():
	parser = argparse.ArgumentParser(description='Filter ace files based on a fasta file')
	parser.add_argument('-f', dest="input_fasta", nargs="*", type=str, help="...")
	parser.add_argument('-t', dest="taxsbp_bins", type=str, help="...")
	parser.add_argument('-g', dest="gzip", default=False, action='store_true', help="...")
	parser.add_argument('-o', dest="output_folder", type=str, help="...")
	args = parser.parse_args()


	taxsbp_write_files(args.input_fasta, args.taxsbp_bins, args.gzip, args.output_folder)
	
def taxsbp_write_files(files, taxsbp_bins_file, gzip_output, output_folder):

	### TODO - double check appending on files, remove them first if possible (specially tmp) ###
	bins = {}
	with open(taxsbp_bins_file, 'r') as file:
		for line in file:
			seqid, _, _, binid = line.rstrip().split("\t")
			bins[seqid] = binid

	for file in files:
		handler = gzip.open(file, 'rt') if file.endswith(".gz") else open(file,'r')
		for header, sequence in fasta_read_generator(handler):
			if header.split(' ', 1)[0] in bins:
				outfilename = output_folder + bins[header.split(' ', 1)[0]] + ".fna"
				outfile = gzip.open(outfilename+".gz", 'at') if gzip_output else open(outfilename, 'a')
				outfile.write('>' + header + '\n' + sequence)
				outfile.close()

def fasta_read_generator(file_handler):
    """
    This function returns a generator that yields name (without the leading >) and sequence tuples from a fasta file of
    which the handler is given.
    
    :param file_handler: file object wrapper for the input fasta file 
    :return: generator yielding name, sequence tuples 
    """
    seq = []
    name = ''
    for line in file_handler:
        if line[0] == '>':
            sequence = ''.join(seq)
            if name:  # only yield when we already have all data for the first sequence
                yield name, sequence
            name = line.rstrip()[1:]  # omitting the leading >
            seq = []
        else:
            seq += [line]#.rstrip()] # keep line breaks
    sequence = ''.join(seq)
    yield name, sequence  # don't forget the last sequence

if __name__ == '__main__':
	main()
