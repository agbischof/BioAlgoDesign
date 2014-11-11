#
def PatternCount(text, pattern):
	"""Count the number of occurances of k-mer, Pattern, in string text"""
	count = 0
	for i in range(len(text)-len(pattern)+1):

		if text[i:i+len(pattern)] == pattern:
			count += 1
	return count


def FrequentWords(text, k):
	"""Find the most frequent k-mer of length k in string text"""
	freqPatt = []
	count = []
	for i in range(len(text)-k+1):
		pattern = text[i:i+k]
		count.append(PatternCount(text, pattern))

	maxCount = max(count)

	for i in range(len(text)-k+1):
		if count[i] == maxCount:
			if text[i:i+k] not in freqPatt:
				freqPatt.append(text[i:i+k])

	return sorted(freqPatt) 


def PatternMatch(text, pattern):
	"""Record all positions in the string text where kmer pattern is found"""
	positions = []

	for i in range(len(text)-len(pattern)+1):
		if text[i:i+len(pattern)] == pattern:
			positions.append(i)

	return positions


def ClumpFinder(genome, k, L, t):
	"""Find distinct kmers of length k which show up t times 
	in window of genome of length L"""
	
	FreqWords = []	
	kmer_dict = {}

	#create a dictionary of all the kmers in genome 
	#the value associated with each kmer is a list
	#the first entry keeps track of the total counts of the kmer
	#the remainder of the list is the index of each occurance
	for i in range(len(genome)-k+1):
		kmer = genome[i:i+k]
		if kmer not in kmer_dict:
			kmer_dict[kmer] = [1,i]
		else:
			kmer_dict[kmer][0] += 1
			kmer_dict[kmer].append(i)

	for kmer in kmer_dict:
		if kmer_dict[kmer][0] >= t and kmer not in FreqWords:
			for n in range(len(kmer_dict[kmer])-t):
				if kmer_dict[kmer][n+t-1] - kmer_dict[kmer][n] <= L:
					if kmer not in FreqWords:
						FreqWords.append(kmer)

	return sorted(FreqWords)
#text = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
#f = open("E-coli.txt", 'r')
#text = f.read()
#kmer = 9
#window = 500
#repeats = 3

text = "CCGACAGGCTAGTCTATAATCCTGAGGCGTTACCCCAATACCGTTTACCGTGGGATTTGCTACTACAACTCCTGAGCGCTACATGTACGAAACCATGTTATGTAT"
kmer = 4
window = 30
repeats = 3

print ClumpFinder(text, kmer, window, repeats)
#f.close()