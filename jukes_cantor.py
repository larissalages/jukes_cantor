import numpy as np
import csv

"""As the matrix is sparse, storing in a square matrix would waste unnecessary memory space,
   so I created this class to store only the necessary distances.
   The class has four methods:
   -add: Stores distance of two sequences
   -getValue: Returns the value of the distance between two sequences
   -printMatrix: Prints at distances from sequences on the console
   -save_to_csv: Saves the information stored in the object in a csv
"""
class distance_matrix:

	def __init__(self):
		global dict1
		dict1 = {}

	def add(self,label_seq1,label_seq2,distance):

		if (label_seq2 not in dict1) or (label_seq2 in dict1 and label_seq1 not in dict1[label_seq2]): #If dict1[label_seq2][label_seq1] existed, it would not be necessary to store the value again
				if label_seq1 not in dict1:
					dict2 = {}
					dict1[label_seq1] = dict2
				dict1[label_seq1][label_seq2] = distance

	def getValue(self,label_seq1,label_seq2):
		if label_seq1 in dict1:
			if label_seq2 in dict1[label_seq1]:
				return dict1[label_seq1][label_seq2]

	def printMatrix(self):
		for key1 in dict1:
			for key2 in dict1[key1]:
				print("Distance %s  &  %s: %.6f\n"%(key1,key2,dict1[key1][key2]))

	def save_to_csv(self,output_file_name):
		with open(output_file_name, 'w') as csvfile:
			fieldnames = ['Sequence 1', 'Sequence 2', 'Distance']
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

			writer.writeheader()
			for key1 in dict1:
				for key2 in dict1[key1]:
					writer.writerow({'Sequence 1': key1, 'Sequence 2': key2,'Distance': dict1[key1][key2]})

			csvfile.close()


#Read input file
def read_file(file_name):

	list_labels = []
	list_seq = []

	try:
		with open(file_name, 'r') as f:       

			count = 0

			for line in f:
				if (count % 2 == 0):
					list_labels.append(line)
				else:
					list_seq.append(line)
				count+=1

			f.close()
		
	except IOError: 
		print "Could not read file"

	return list_labels, list_seq   

#Remove the characters ">" and "\n"
def remove_inval_char(label_seq1,label_seq2):
	label_seq1 = label_seq1.replace(">", "")
	label_seq1 = label_seq1.replace("\n", "")

	label_seq2 = label_seq2.replace(">", "")
	label_seq2 = label_seq2.replace("\n", "")

	return label_seq1, label_seq2


#Calculate Jukes Cantor distance
def calculate_jukes_cantor(seq1,seq2):

	length = 0 #length counter
	diff = 0	#difference counter
	nucleotides = 'ACGT'

	#find percentage of characters that are different between two sequences
	for i in range(min(len(seq1), len(seq2))):
		if (seq1[i] in nucleotides) and (seq2[i] in nucleotides):
			length += 1
			if seq1[i] != seq2[i]:
				diff += 1

	#The equation does not exist for values of p greater than 3/4, neither for values of length == 0
	if(length != 0): 
		p = float(diff)/length
		if(p < 3.0/4):
			d = -3.0/4*np.log(1 - 4.0/3*p)
			if d == -0.0:
				d = d*(-1)
		else:
			d = np.inf
	else:
		d = np.inf

	return d


# Calculate Jukes Cantor distance for each pair of sequences
def calculate_matrix(list_seq,list_labels):
	dist_matrix = distance_matrix()

	for i in range(len(list_seq)):
		for j in range(len(list_seq)):
			if i != j:
				dist = calculate_jukes_cantor(list_seq[i],list_seq[j])
				if dist < 0.04: #discard values greater than 4%
					label1,label2 = remove_inval_char(list_labels[i],list_labels[j])
					dist_matrix.add(label1,label2,dist)

	return dist_matrix


def main():

	#Read input file
    list_labels, list_seq = read_file("Public_Eacles_Morpho_Aglia.fas")

    #Calculate the distance between these sequences
    dist_matrix = calculate_matrix(list_seq,list_labels)

    #Save the output in a CSV file
    dist_matrix.save_to_csv('output_matrix.csv')


if __name__ == "__main__":
    main()