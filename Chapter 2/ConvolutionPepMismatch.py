

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

def LinearPepSpect(peptide):
	"""output mass spectrum of linear peptide, peptide given in numbers"""
	spectrum = [0]
	PrefixMass = [0]

	for i in range(len(peptide)):
		PrefixMass.append(peptide[i] + PrefixMass[i])

	for i in range(len(peptide)):
		for j in range(i+1, len(peptide)+1):
			spectrum.append(PrefixMass[j] - PrefixMass[i])

	return sorted(spectrum)

def CycloPepScore(peptide, spectrum, aaMass):
	"""Return the score of peptide against spectrum (# of masses shared)"""
	TheoSpect = CyclicPepSpect(peptide, aaMass)
	score = 0

	spect = list(spectrum)
	for mass in TheoSpect:
		if mass in spect:
			score += 1
			spect.remove(mass)

	return score

def LinearPepScore(peptide, spectrum):
	"""Return the score of peptide against spectrum (# of masses shared)"""
	TheoSpect = LinearPepSpect(peptide)
	score = 0

	spect = list(spectrum)
	for mass in TheoSpect:
		if mass in spect:
			score += 1
			spect.remove(mass)

	return score

def expand(peptides, aaMasses):
	"""Expand all peptides in list by 1 amino acid, making each possible combination"""
	newpeptides = []
	for pep in peptides:
		for mass in aaMasses:
			if pep == 0:
				newpeptides.append([pep] + [mass])
			else:
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

def Trim(Board, spectrum, n):
	Scoreboard = []
	KeepPeps = []

	for pep in Board:
		Scoreboard.append(LinearPepScore(pep, spectrum))

	#print Scoreboard
	Scoreboard.sort(reverse = True)
	#print Scoreboard
	if len(Scoreboard) >= n:
		scoreCut = Scoreboard[n-1]
		for pep in Board:
			if LinearPepScore(pep, spectrum) >= scoreCut:
				KeepPeps.append(pep)
	else:
		KeepPeps = Board

	return KeepPeps

def LeaderboardCycSeq(spectrum, n, aaMass):
	Leaderboard = expand([0], aaMass)

	LeaderPeptide = ""
	removePeps = []

	while len(Leaderboard) > 0:
		#print Leaderboard
		removePeps = []
		for peptide in Leaderboard:
			print 'pep:', peptide
			#print 'spec:', spectrum
			if sum(peptide) == max(spectrum):
				if LinearPepScore(peptide, spectrum) > LinearPepScore(LeaderPeptide, spectrum):
					LeaderPeptide = peptide
			elif sum(peptide) > max(spectrum):
				removePeps.append(peptide)

		for pep in removePeps:
			Leaderboard.remove(pep)
		Leaderboard = Trim(Leaderboard, spectrum, n)
		if len(Leaderboard) > 0:
			Leaderboard = expand(Leaderboard, aaMass)

	return LeaderPeptide		


def SpectralConv(spectrum):
	"""Find the convolution of spectrum"""
	convol = []

	for i in range(len(spectrum)):
		for j in range(i+1, len(spectrum)):
			value = spectrum[j] - spectrum[i]
			if value > 0:
				convol.append(spectrum[j] - spectrum[i])
	return sorted(convol)

def ConvolutionCyclopeptideSeq(M, N, spectrum):
	"""determine the cyclic peptide which is the best fit for the spectrum made up of the top M elemnets 
	based on the convolution of spectrum that fall between 57 and 200, where leaderboard is restricted to N"""

	aaDict = {}
	aaMass = []
	spectrum.sort()
	Convolution = SpectralConv(spectrum)
	for mass in Convolution:
		if mass >= 57 and mass <= 200:
			if mass in aaDict:
				aaDict[mass] += 1
			else:
				aaDict[mass] = 1

	if len(aaDict) >= M:			
		KeepFreqCut = sorted([aaDict[k] for k in aaDict])[M]
		AAs = aaDict.keys()
		for element in AAs:
			if aaDict[element] >= KeepFreqCut:
				aaMass.append(element)
	else:
		aaMass = aaDict.keys()
	#print 'Mass: ', aaMass
	peptides = LeaderboardCycSeq(spectrum, N, aaMass)

	return peptides




aaMass = {}

with open("aa_integer_mass_table.txt") as f:
	for line in f:
		#print line
		(key, val) = line.split(" ")
		val = val.rstrip('\n')
		aaMass[key] = val

with open("TestSpect.txt") as f:
	for line in f:
		spectrum = line.split()
		if spectrum:
			spectrum = [int(i) for i in spectrum]

with open("PepList.txt") as f:
	for line in f:
		Peps = line.split(' ')
		if Peps:
			Peps = [str(i) for i in Peps]



spectrum.sort()
print ConvolutionCyclopeptideSeq(20, 373, spectrum)

#print expand([[128], [121], [150]], [91, 93, 99])

