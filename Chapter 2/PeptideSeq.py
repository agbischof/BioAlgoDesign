
def CycloSeqNum(l):
	"""Calculate the number of possible peptides of a cyclic peptide of length, l"""
	return l*(l-1)


def LinearPepSpect(peptide, aa_Mass):
	"""output mass spectrum of linear peptide, peptide"""
	spectrum = [0]
	PrefixMass = [0]

	for i in range(len(peptide)):
		PrefixMass.append(int(aa_Mass[peptide[i]]) + PrefixMass[i])

	for i in range(len(peptide)):
		for j in range(i+1, len(peptide)+1):
			spectrum.append(PrefixMass[j] - PrefixMass[i])

	return sorted(spectrum)

def CyclicPepSpect(peptide, aa_Mass):
	"""output mass spectrum of cyclic peptide, peptide"""
	spectrum = [0]
	PrefixMass = [0]

	for i in range(len(peptide)):
		PrefixMass.append(int(aa_Mass[peptide[i]]) + PrefixMass[i])
	PrefixMass = PrefixMass[1:]
	peptideMass = PrefixMass[-1]
	for i in range(len(peptide)-1):
		for j in range(i+1, len(peptide)):
			spectrum.append(PrefixMass[j] - PrefixMass[i])
			spectrum.append(peptideMass - (PrefixMass[j] - PrefixMass[i]))
	spectrum.append(peptideMass)
	return sorted(spectrum)


def CountPeptides(m):
	pass

def CountSubpeptides(n):
	"""Output the number of subpeptides of a linear peptide of length n"""
	count = n*(n-1)/2.0

	return count

def expand(peptides):
	"""Expand all peptides in list by 1 amino acid, making each possible combination"""
	aaMasses = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
	newpeptides = []
	for pep in peptides:
		for mass in aaMasses:
			newpeptides.append(pep + [mass])
	return newpeptides

def MassLinear(peptide, TheSpectrum):
	"""Give mass spectrum for given peptide listed in numbers representing the mass of individual amino acids
	return True if spectrum of peptide is consistent with TheSpectrum, False if not"""
	PrefixMass = [0]

	for i in range(len(peptide)):
		PrefixMass.append(peptide[i] + PrefixMass[i])
	#print PrefixMass
	for i in range(len(peptide)+1):
		for j in range(i+1, len(peptide)+1):
			if PrefixMass[j] - PrefixMass[i] in TheSpectrum:
				TheSpectrum.remove(PrefixMass[j] - PrefixMass[i])
			else:
				return False
	return True

def MassCyclic(peptide):
	"""output mass spectrum of cyclic peptide, peptide (listed in numbers representing indiv amino acids"""
	spectrum = [0]
	PrefixMass = [0]

	for i in range(len(peptide)):
		PrefixMass.append(peptide[i] + PrefixMass[i])
	PrefixMass = PrefixMass[1:]
	peptideMass = PrefixMass[-1]
	for i in range(len(peptide)-1):
		for j in range(i+1, len(peptide)):
			spectrum.append(PrefixMass[j] - PrefixMass[i])
			spectrum.append(peptideMass - (PrefixMass[j] - PrefixMass[i]))
	spectrum.append(peptideMass)
	return sorted(spectrum)


def CyclopeptideSeq(spectrum):
	FinalSpect = []
	spectrum.sort
	peptides = [[57], [71], [87], [97], [99], [101], [103], [113], [114], [115], [128], [129], [131], [137], [147], [156], [163], [186]]

	while len(peptides) > 0:
		#print 'first in while:', peptides
		#print peptides
		removePep = []
		for pep in peptides:
			
			#print pep
			#print MassLinear(pep, list(spectrum))
			
			if MassLinear(pep, list(spectrum)):
				if MassCyclic(pep) == spectrum:
					FinalSpect.append(pep)
					removePep.append(pep)
					#print 'final match:', FinalSpect
			else:
				removePep.append(pep)
		for pep in removePep:
			peptides.remove(pep)

		peptides = expand(peptides)

	return FinalSpect





aaMass = {}

with open("aa_integer_mass_table.txt") as f:
	for line in f:
		#print line
		(key, val) = line.split(" ")
		val = val.rstrip('\n')
		aaMass[key] = val

#print ' '.join(map(str, CyclicPepSpect('GMWQGLTNHERIK', aaMass)))
#print CountSubpeptides(17059)
#spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
with open("TestSpect.txt") as f:
	for line in f:
		spectrum = line.split()
		if spectrum:
			spectrum = [int(i) for i in spectrum]

#print spectrum
Seqs = CyclopeptideSeq(spectrum)
SeqFinal = []

for pep in Seqs:
	SeqFinal.append('-'.join(str(aaMass) for aaMass in pep))

print ' '.join(SeqFinal) 

#print ' '.join(str(pep) for pep in CyclopeptideSeq(spectrum))
#print MassLinear([97], list(spectrum))
#print MassLinear([113], list(spectrum))
#print spectrum
#print MassLinear([113], list(spectrum))
