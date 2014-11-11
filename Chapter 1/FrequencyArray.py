#Use frequency arrays to find patterns in DNA 
#ClumpFinder finds kmers (length k) that repeat at least t times in a window of DNA (length L)

def SymbolToNumber(symbol):
	"""Convert base into alphabetically ordered number"""
	baseConvert = {'A':0, 'C':1, 'G':2, 'T':3}
	number = baseConvert[symbol]
	return number

def NumberToSymbol(number):
	"""Convert base into alphabetically ordered number"""
	baseConvert = {0:'A', 1:'C', 2:'G', 3:'T'}
	symbol = baseConvert[number]
	return symbol


def NumberToPattern(index, k):
	"""Convert number to pattern of kmer length k given that the index 
	is a number ordered of all possible combinations of kmers length k ordered alphabetically"""
	pattern = []

	if k == 1:
		return NumberToSymbol(index)
	else:
		PrefixPattern = NumberToPattern(int(index/4), k-1)
		symbol = NumberToSymbol(index % 4)
		return PrefixPattern + symbol


def PatternToNumber(Pattern, debug = False):
	"""Convert pattern into indexed number of all patterns of the same size ordered alphabetically"""
	number = 0
	if len(Pattern) == 1:
		return SymbolToNumber(Pattern)
	else:
		number = number + 4 * PatternToNumber(Pattern[:-1]) + SymbolToNumber(Pattern[-1])  
		if debug:	
			print Pattern[:-1]
			print Pattern[-1]
			print number
			print
		return number

def ComputingFrequencies(text, k, debug = False):
	"""Computes frequency array containing counts of the total occurances of kmers (length k) 
	in a string (text). The numbers are ordered as if the kmers were ordered alphabetically"""

	kmer_dict = {}
	FreqArray = []
	for i in range(4**k):
		FreqArray.append(0)
	if debug:
		print "freqArray is length:", len(FreqArray)
		
	for i in range(len(text)-k+1):
		Pattern = text[i:i+k]
		j = PatternToNumber(Pattern)
		FreqArray[j] = FreqArray[j] + 1
		if debug:
			print j
			print Pattern
	return FreqArray	


def ClumpFinder(Genome, k, L, t):
	"""Finds kmers (length k) that repeat at least t times in a window of DNA (length L)""" 

	FreqPatterns = []
	Clump = []
	
	for i in range(4**k):
		Clump.append(0)
	
	Text = Genome[0:L]
	print len(Text)
	FreqArray = ComputingFrequencies(Text, k)
	
	for i in range(4**k):
		if FreqArray[i] >= t:
			Clump[i] = 1

	for i in range(1, len(Genome)-L):
		FirstPattern = Genome[i-1:i+k-1]
		j = PatternToNumber(FirstPattern)
		FreqArray[j] = FreqArray[j] - 1

		LastPattern = Genome[i+L-k:i+L]
		j = PatternToNumber(LastPattern)
		FreqArray[j] = FreqArray[j] + 1

		if FreqArray[j] >= t:
			Clump[j] = 1

	for i in range(4**k):
		if Clump[i] == 1:
			Pattern = NumberToPattern(i,k)
			FreqPatterns.append(Pattern)
	return FreqPatterns

if __name__ == "__main__":
	"""Test Code"""
	# print PatternToNumber('ATGCAA')
	# print NumberToPattern(5437, 8)

	#print ' '.join(str(x) for x in ComputingFrequencies("ATC", 8))
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
