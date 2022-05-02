import argparse

def replace_str_index(text,index=0,replacement=''):
    return '%s%s%s'%(text[:index],replacement,text[index+1:])

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('unwrap_fasta', help='')
	args = parser.parse_args()

	# Grabs name of file.
	fasta_name = str(args.unwrap_fasta) + "_taylor.fasta"
	new_fasta = open(fasta_name, "w")

	# Grabs only the fasta we want (not the reference fasta),
	# with the bases 201-29740 (per sarscov2_masterfile.txt).
	line_num = 0
	n_indexes = []
	for line in open("repositioned_unwrap.fasta"):
		# Finds correctly masked positions in reference
		if(line_num == 1):
			n_indexes = find(line, "n")
		# Writes fasta header with real name
		elif(line_num == 2):
			fasta_header = ">" + args.unwrap_fasta + "\n"
			new_fasta.write(fasta_header)
		# Masks positions by ref masked
		elif(line_num==3):
			line_masked = line
			for index in n_indexes:
				if(line_masked[index]!="-"):
					line_masked = replace_str_index(line_masked,index,"n")
			line_trimmed = line_masked[201:29741]
			final_line = line_trimmed.replace("-","")
			new_fasta.write(final_line)
		line_num = line_num + 1