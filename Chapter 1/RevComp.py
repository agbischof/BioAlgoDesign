# Find the reverse complement of a DNA string

def RevComp(pattern):
	basecomplement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	letters = list(pattern)
	letters = [basecomplement[base] for base in letters]

	return ''.join(letters)[::-1]

#print RevComp("AAAACCCGGT")


if __name__ == "__main__":
	sequence = "TTGTGTC"

	print RevComp(sequence)

