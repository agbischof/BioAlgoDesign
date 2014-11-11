
import random

#random string modeling
def ProbExact(l, Letters, repeats):
	"""Determine the probability that 2 strings of length, l, made of of letters from list Letters are identical"""
	
	equal = 0
	for rep in range(repeats):
		str1 = ""
		str2 = ""

		for i in range(l):
			str1 = str1 + Letters[random.randrange(0,4, 1)]
			str2 = str2 + Letters[random.randrange(0,4,1)] 

		if str1 == str2:
			equal += 1

	print equal, '/', repeats
	return float(equal)/repeats


print ProbExact(3, ['A', 'C', 'T', 'G'], 100000)

print 2/(float(4)**9*2)
print (float(4)**9*2)


def Pr(N, A, pattern, t):

	include = 0
	array_text = []

	for  