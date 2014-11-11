#Comput Hamming distance between 2 strings

def HammingDist(str1, str2):
	"""Compute the Hamming Distance between the 2 strings, str1 and str2
	Hamming Distance is number of mismatches between the strings"""
	Hdist = 0
	for i, base in enumerate(str1):
		if base != str2[i]:
			Hdist += 1

	return Hdist
#test code
#text1 = "ATGGAGCACCCGCAATGCATGGTAGAAATTTGCACCGTCAACACAAGGACGCCCGCTCTCCAATGATCTTCGTAGCTGATTCAATTTGTACGGTCTCTATCTGCCACGGATCGAGAGCGTGTGCAGGGTATCTAGGTGCTTGGTCATCCTACTAAGATACACTAAGACCTTGCAGCAAAAACCACGGAGGCATCGTCTTTTATGTGGTGTTATCCTGCATCATGAGGCCAAAGCATCGCTCAGTTTCTCTAGGCGCGCGGTCGCGTCCCGCCAATTATGACCCGGTCCGTTTTGCGTTGGACTCGAAAACTCCTGAGGCCATCTGAACGGCCATCGTAAAGTTCAATAACACACGCGCTCGCCAACCTGGACCCATGTTAATCACTCGGCCCACACGTATGTGCCTCCCCTCGTCACAATGCACCTTACGTTTGCGTCGCTCTTTCGTTGCTCAAGGCAGGCACTATGAGGGTAGTCGGGTTGCGAGCGTGCCTCACCCCTGAGTCTCCTGCCCTTCGTACTCGACGGATACTGTGGTGAATTGCCAACTCAATTACAATGATCTGTTGTGTTCTAATGTCTGCCACCGGTACGTCCTTGCTGTTGGTATAAGTCTCCGGTTTACCTGCAGGATAGACGAAGTATGTACCGACAAAGAAGGACGGACTCTCTAGGACGCTTTCATTGAGACCTGTTTCGCATGCTTTGGCCGCAGACCAAAGACGCCCTCGAGTCTCAGCCGGTGCAGACGAAAAAACCGGCCCCTGGAGATGCTTCTGTCGTCAAAACTAGATGAGTAGTGGTGCAAGGAGAGTCAACCCGGTAATTACGAAGAACTTTAGTCTCGTGTGTCGGATCCGCCCCAAGTTCTCTCGCACCGGTTAAGTCATTGTCGGAGCCATGGCTGAGGGGCCCGAGGCTTCTTTGCTTCAATACGCTACCACTAACTCGTCGACTCTGATTGGCTAGGGAATCCAGCGCCACGTACAGAGTAAGGTACCGTTCAGGGAGACGGTAATAGCACGCCCATCGCGCCGCTCATTAAATTCAGCGCATACGCCAAGACCAGACGCAATTGTTAGTGTACTCATTCGCTGACCGGTAGATTACAACGTATGCGTAGAGCCTGCTAAT"
#text2 = "CTCCGACCAAACCTGTAGTGACCAAGCCTAAGCCACCTAGTGTTATAAGCACCGATGGCTAGCAAGTGACACCGTGGGGCTACATCTCACAAACGCTCTGCCGCGAACCATCTAACAACCAGGGTCCCCGAGGCCAGGGTTGTGATTCATATGCGAACGGAGTAGGATTCCTCCAATGGGTAGGAAGAGCCGTTGACCAGGCTTGTGTAGAGCGTCATAACTTTAATGTTATACGGTGCTCAGACCAGAAGGCTAGTGGCCTGGTCTATGTGAAATTGCGCTACTTCACTCTCGGGGGCATTTCGGCCGTCCCTCGAGGGATAGTGTAGCAAGGGTACGCATACGCAAATTGATCCCGTTGCCGGTTTGCTGCCCTGACGTTACTGCGTGCTTACGTTATTCCAAAGTTGGATTGCGTACCGACCGGGACAGGAATGATTCAGAGCGCAAGCTTCAGACGGCGGATTAGATAATCTTTATTGTGCATACTGCTAATTGGATGGTGTATGTTCATATAAGCGGCCTCGTCCTGATCCCCCGGTCTCGACGAGGGGCCCTGCGGGCACTGAGATGTATAGAGTGACTGGTACGCCTGCTGCAAGCCTAGCGCCCCGTGAGAATAATGTTCCCATATAACGCATTTAAGGAGCGAAGATCATCGCTTAAATATTCTTTTAATGCTCACAACTGCGTCGGCTGGACTCGCGGATGTAGTACCGTTGGGCTTCAGGCCTAATCACGACTGCGCGTAACTCGAATTGGGCAGCAAGATCTCCTTGGGTTCGTCTGGAATACTTGCCAACAAGCATCCTATTCCAATACTTCGCTTCTGAAATACCGTAGACCCTGGCTGGATTTTCTTGGCGCCGCTTCGGTTCAGAGCTGCACTTGGCCACACAATAAGCACTCTGCCAATCCGGGAAACGTATGCCGCCTCCGGCTTAGGTCGCAACACTACAACGTTGGTGGACTAGTCAACCCACTTTTTCTATGTTCTGCCCGGACATAAAGCGTTAAAGAGACTCAACTCCTTGAACCATTCAAGATACGTGTTATTACTAGTTCCATTCGTCCTGCACGCTACGCCTTAGTTCGGATGACCCCAAGTTAACGCACATTTTGCGTGTCTATCAT"

#print HammingDist(text1, text2)




def ApproxStrMatch(pattern, text, d, debug = False):
	"""Find all approximate occurrences of pattern in string with at most d mismatches"""
	count = 0
	index = []
	for i in range(len(text)-len(pattern)+1):
		if debug:
			print text[i:i+len(pattern)]
			print HammingDist(text[i:i+len(pattern)], pattern)
		if HammingDist(text[i:i+len(pattern)], pattern) <= d:
			count += 1
			index.append(i)
	return index

def CountAppStrMatch(pattern, text, d, debug = False):
	"""Find all approximate occurrences of pattern in string with at most d mismatches"""
	count = 0
	if debug:
		print len(text)-len(pattern)+1
	for i in range(len(text)-len(pattern)+1):
		if debug:
			print text[i:i+len(pattern)]
			print HammingDist(text[i:i+len(pattern)], pattern)
		if HammingDist(text[i:i+len(pattern)], pattern) <= d:
			count += 1
	return count

def FreqWordsMisMatch(text, k, d):
	AllWords = []	
	kmer_dict = {}
	AllMisMatch = []

	#create a dictionary of all the kmers in genome 
	#the value associated with each kmer is a list
	#the first entry keeps track of the total counts of the kmer
	#the remainder of the list is the index of each occurance
	for i in range(len(text)-k+1):
		kmer = text[i:i+k]
		if kmer not in kmer_dict:
			kmer_dict[kmer] = [1,i]
			AllWords.append(kmer)
		else:
			kmer_dict[kmer][0] += 1
			kmer_dict[kmer].append(i)

	for kmer in AllWords:
		for pattern in AllWords:
			if HammingDist(kmer, pattern) <=d and kmer not in AllMisMatch:
				AllMisMatch.append(kmer)

	return AllMisMatch

#Calculate all mutations within the specified distance of word
#figure out how this works or design your own..
def mutations(word, hamming_distance): 
	for indices in itertools.combinations(range(len(word)), hamming_distance): 
		for replacements in itertools.product('ATCG', repeat=hamming_distance): 
			mutation = list(word) 

		for index, replacement in zip(indices, replacements): 
			mutation[index] = replacement 
		yield "".join(mutation)

if __name__ == "__main__":
	#kmer = "TCCTG"
	genome = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
	mismatch = 1
	k = 4

	#print FreqWordsMisMatch(genome, k, mismatch)
	#print CountAppStrMatch(kmer, genome, mismatch)
	#print ' '.join(str(x) for x in ApproxStrMatch(kmer, genome, mismatch))

	text1 = "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
	text2 = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
	print "HammingDistance:", HammingDist(text1, text2)

	#pattern = "TGT"
	#text = "CGTGACAGTGTATGGGCATCTTT"
	#d = 1

	#print CountAppStrMatch(pattern, text, d)


