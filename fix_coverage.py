import os
import sys
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('sample_name', help='sample name')
	args = parser.parse_args()

	# Grabs name of file.
	summmary_csv = str(args.sample_name) + "_summary.csv"
	summary_csv_fixed = open(str(args.sample_name) + "_summary_fixed.csv","w+")

	line_num = 0
	n_percent = 100
	for line in open(str(args.sample_name) + "_taylor.fasta"):
		if(line_num==1):
			sequence = line
			counter = sequence.lower().count('n')
			length = len(sequence)
			n_percent = (counter/length*100)
		line_num +=1

	line_num = 0
	for line in open(summmary_csv):
		if(line_num==0):
			summary_csv_fixed.write(line)
		if(line_num==1):
			sample_data = line.split(",")
			new_sample_data = sample_data
			new_sample_data[-1] = n_percent

			new_sample_line = ','.join(map(str,new_sample_data))
			summary_csv_fixed.write(new_sample_line)
		line_num+=1