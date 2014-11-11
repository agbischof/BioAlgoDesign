

def HammingDist(str1, str2):
	"""Compute the Hamming Distance between the 2 strings, str1 and str2
	Hamming Distance is number of mismatches between the strings"""
	Hdist = 0
	for i, base in enumerate(str1):
		if base != str2[i]:
			Hdist += 1

	return Hdist

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

def ImmediateNeighbors(pattern):
	""""""
	Neighborhood = [pattern]

	for i in range(len(pattern)):
		symbol = pattern[i]
		for x in 'ATCG':
			if x != symbol:
				Neighbor = pattern[:i] + x + pattern[i+1:]
				Neighborhood.append(Neighbor)

	return Neighborhood


def Neighbors(pattern, d):
	"""Create a list of all strings within d Hamming distance from pattern"""
	if d == 0:
		return pattern
	if len(pattern) == 1:
		return {'A', 'C', 'G', 'T'}

	Neighborhood = []
	SuffixNeighbors = Neighbors(pattern[1:], d)
	for text in SuffixNeighbors:
		if HammingDist(pattern[1:], text) < d:
			for nucleotide in 'ATCG':
				Neighborhood.append(nucleotide+text)
		else:
			Neighborhood.append(pattern[0]+text)

	return Neighborhood


def CountAppStrMatch(pattern, text, d, debug = False):
	"""Find all approximate occurrences of pattern in string with at most d mismatches"""
	count = 0
	for i in range(len(text)-len(pattern)+1):
		if debug:
			print text[i:i+len(pattern)]
			print HammingDist(text[i:i+len(pattern)], pattern)
		if HammingDist(text[i:i+len(pattern)], pattern) <= d:
			count += 1
	return count

def RevComp(pattern):
	basecomplement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	letters = list(pattern)
	letters = [basecomplement[base] for base in letters]

	return ''.join(letters)[::-1]

#print RevComp("AAAACCCGGT")

def FrequentWordsMismatches(text, k, d):
	"""Find the most frequent kmers with mismatches up to d in text"""
	FrequentPatterns = []
	Close = []
	FrequencyArray = []

	for i in range(4**k-1):
		Close.append(0)
		FrequencyArray.append(0)

	for i in range(len(text)-k+1):
		Neighborhood = Neighbors(text[i:i+k], d)
		for pattern in Neighborhood:
			index = PatternToNumber(pattern)
			Close[index] = 1

	for i in range(4**k - 1):
		if Close[i] == 1:
			pattern = NumberToPattern(i,k)
			FrequencyArray[i] = CountAppStrMatch(pattern, text, d)

	maxCount = max(FrequencyArray)
	print maxCount

	for i in range(4**k - 1):
		if FrequencyArray[i] == maxCount:
			pattern = NumberToPattern(i,k)
			FrequentPatterns.append(pattern)

	return FrequentPatterns 

def FindingFreqWordsBySort(text, k, d, debug = False):
	"""Find the most frequent kmers with mismatches up to d in text"""
	FrequentPatterns = []
	Neighborhoods = []
	index = []
	count = []

	for i in range(len(text)-k+1):
		Neighborhoods.append(Neighbors(text[i:i+k], d))

	if debug:
		#print 'Neighborhoods:', Neighborhoods
		pass

	NeighborhoodArrays = [item for sublist in Neighborhoods for item in sublist] #[[x] for x in Neighborhoods]
	if debug:
		#print 'NeighborhoodArrays:', NeighborhoodArrays
		print 'Len NeighborhoodArrays:', len(NeighborhoodArrays)

	for i in range(len(NeighborhoodArrays)):
		pattern = NeighborhoodArrays[i]
		if debug:
			print pattern, ' ', PatternToNumber(pattern)
			pass
		index.append(PatternToNumber(pattern))
		count.append(1)
	
	
	SortedIndex = sorted(index)
	if debug:
		print 'sorted index:', SortedIndex
		print 'length of count:', len(count)
		print 'count pre:', count
		pass
	for i in range(len(SortedIndex)-1):
		if SortedIndex[i] == SortedIndex[i+1]:
			count[i+1] = count[i] + 1
	if debug:
		print count
			
	maxCount = max(count)
	print maxCount
	for i in range(len(SortedIndex)):
		if count[i] == maxCount:
			pattern = NumberToPattern(SortedIndex[i], k)
			FrequentPatterns.append(pattern)
	return FrequentPatterns



def FreqWordsMismatchRevCompBySort(text, k, d, debug = False):
	"""Find the most frequent kmers with mismatches up to d in text"""
	FrequentPatterns = []
	Neighborhoods = []
	index = []
	count = []
	pattern2 = RevComp(text)

	for i in range(len(text)-k+1):
		Neighborhoods.append(Neighbors(text[i:i+k], d))

	for i in range(len(pattern2)-k+1):
		Neighborhoods.append(Neighbors(pattern2[i:i+k], d))

	if debug:
		#print 'Neighborhoods:', Neighborhoods
		pass

	NeighborhoodArrays = [item for sublist in Neighborhoods for item in sublist] #[[x] for x in Neighborhoods]
	if debug:
		#print 'NeighborhoodArrays:', NeighborhoodArrays
		print 'Len NeighborhoodArrays:', len(NeighborhoodArrays)

	for i in range(len(NeighborhoodArrays)):
		pattern = NeighborhoodArrays[i]
		if debug:
			print pattern, ' ', PatternToNumber(pattern)
			pass
		index.append(PatternToNumber(pattern))
		count.append(1)
	
	
	SortedIndex = sorted(index)
	if debug:
		print 'sorted index:', SortedIndex
		print 'length of count:', len(count)
		print 'count pre:', count
		pass
	for i in range(len(SortedIndex)-1):
		if SortedIndex[i] == SortedIndex[i+1]:
			count[i+1] = count[i] + 1
	if debug:
		print count
			
	maxCount = max(count)
	print maxCount
	for i in range(len(SortedIndex)):
		if count[i] == maxCount:
			pattern = NumberToPattern(SortedIndex[i], k)
			FrequentPatterns.append(pattern)
	return FrequentPatterns



#str1 = "ACCCGCTGTCCCCAGACCCGCTCGGCCAGGTCCGCTCCAGCGGACCCGTCCGTCCGCTGTCCGCTGCTACCCGCTCGGACCCCCAGGTCCCCAGCGGCGGGCTGCTACCCCGGGCTGCTACCCGTCCACCCGTCCGCTCCAGGCTGCTACCCGTCCACCCACCCGTCCCCAGGTCCGCTGTCCCCAGGCTCCAGCGGCCAGGCTGCTCGGACCCGTCCGCTGTCCGTCCGTCCCGGCCAGGTCCCCAGGTCCGTCCACCCGCTGCTGTCCCCAGACCCCCAGCGGACCCGCTCCAGCCAGCGGCGGACCCGTCCCCAGACCCGCTGTCCGTCCGCTGTCCACCCGTCCACCCGCTACCCCCAGGCTGCTACCCGTCCCCAGGCTCGGACCC"
#k = 9
#d = 3
#answers ['CCCCCGCCC', 'CCCCGGCCC']

#str1 = "CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC"
#k = 10
#d = 2
#answers: ['GCACACAGAC', 'GCGCACACAC']

#str1 = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
#k = 4
#d = 1
#answers ['ATGC', 'ATGT', 'GATG']

#print PatternToNumber(str1)
#print FindingFreqWordsBySort(str1, k, d)

#print FrequentWordsMismatches(str1, k, d)


#mismatch and reverse complement
#str1 = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
#k = 4
#d = 1
#answers: ['ACAT', 'ATGT']


#str1 = "CTTGCCGGCGCCGATTATACGATCGCGGCCGCTTGCCTTCTTTATAATGCATCGGCGCCGCGATCTTGCTATATACGTACGCTTCGCTTGCATCTTGCGCGCATTACGTACTTATCGATTACTTATCTTCGATGCCGGCCGGCATATGCCGCTTTAGCATCGATCGATCGTACTTTACGCGTATAGCCGCTTCGCTTGCCGTACGCGATGCTAGCATATGCTAGCGCTAATTACTTAT"
#k = 9
#d = 3
#answers: ['AGCGCCGCT', 'AGCGGCGCT']


#str1 = "CATTCTCCCATCGCCGCCGCCCCGCCGCAGCATCCCATCGCCCAGTCTCATCATTCTCATCGCTCTTCTTCTCGCCATCATTCTCCCCCGCTCTAGAGCGCTCTCCCATCGCCCCATCCTCTTCTCGCCATCGCCGCCCCGCTCTAGCGCCATCCCGCCATCGCCATCGCCCCATAGCCCCCGCCATCGCAGAGCGCAGAGCATCATTCTAGCATAGCATCGCTCTCCAGTCTTCTCAT"
#k = 8
#d = 2
#answers: CCCCCCCC GGGGGGGG

#print ' '.join(FreqWordsMismatchRevCompBySort(str1, k, d))


#pattern = "TGT"
#text = "CGTGACAGTGTATGGGCATCTTT"
#d = 1

#print CountAppStrMatch(pattern, text, d)
pattern = "CCCC"
d = 3
print Neighbors(pattern, d)
print len(Neighbors(pattern, d))