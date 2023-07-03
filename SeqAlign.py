import numpy as np

class Scoring:
	def __init__(self,gap,match,mismatch):
		self.gap = gap
		self.match = match
		self.mismatch = mismatch

	def misMatch(self,x,y):
		if x != y:
			return self.mismatch
		else:
			return self.match

def getAlignedSequences(x,y,matrix,traceBack):
	xSeq = []
	ySeq = []
	i = len(x)
	j = len(y)
	while(i > 0 or j > 0):
		if traceBack[i][j] == 'diag':
			xSeq.append(x[i-1])
			ySeq.append(y[j-1])
			i = i-1
			j = j-1
		elif traceBack[i][j] == 'hor':
			xSeq.append('-')
			ySeq.append(y[j-1])
			j = j-1
		elif traceBack[i][j] == 'ver':
			xSeq.append(x[i-1])
			ySeq.append('-')
			i = i-1
		elif traceBack[i][j] == 'end':
			break
	return xSeq,ySeq

def getMatrix(sizeX,sizeY,gap):
	matrix = []
	for i in range(len(sizeX)+1):
		subMatrix = []
		for j in range(len(sizeY)+1):
			subMatrix.append(0)
		matrix.append(subMatrix)

	for j in range(1,len(sizeY)+1):
		matrix[0][j] = j*gap
	for i in range(1,len(sizeX)+1):
		matrix[i][0] = i*gap
	return matrix

def TraceBack(sizeX,sizeY):
	matrix = []
	for i in range(len(sizeX)+1):
		subMatrix = []
		for j in range(len(sizeY)+1):
			subMatrix.append('0')
		matrix.append(subMatrix)

	for j in range(1,len(sizeY)+1):
		matrix[0][j] = 'hor'
	for i in range(1,len(sizeX)+1):
		matrix[i][0] = 'ver'
	matrix[0][0] = 'end'
	return matrix

###### GLOBAL ALIGNMENT ######
def globalAlign(x,y,score):
	matrix = getMatrix(x,y,score.gap)
	traceBack = TraceBack(x,y)

	for i in range(1,len(x)+1):
		for j in range(1,len(y)+1):
			hor = matrix[i][j-1] + score.gap
			ver = matrix[i-1][j] + score.gap
			diag = matrix[i-1][j-1] + score.misMatch(x[i-1],y[j-1])
			matrix[i][j] = max(hor,ver,diag)
			if matrix[i][j] == hor:
				traceBack[i][j] = 'hor'
			elif matrix[i][j] == ver:
				traceBack[i][j] = 'ver'
			else:
				traceBack[i][j] = 'diag'
	return matrix,traceBack

###### SEMI - GLOBAL ALIGNMENT #####
def semiGlobalAlign(x,y,score):
	matrix = getMatrix(x, y, 0)
	traceBack = TraceBack(x, y)

	for i in range(1, len(x) + 1):
		for j in range(1, len(y) + 1):
			hor = matrix[i][j - 1] + score.gap
			ver = matrix[i - 1][j] + score.gap
			diag = matrix[i - 1][j - 1] + score.misMatch(x[i - 1], y[j - 1])
			matrix[i][j] = max(hor, ver, diag)
			if matrix[i][j] == hor:
				traceBack[i][j] = 'hor'
			elif matrix[i][j] == ver:
				traceBack[i][j] = 'ver'
			else:
				traceBack[i][j] = 'diag'
	return matrix, traceBack

###### LOCAL ALIGNMENT ######
def localAlign(x,y,score):
	matrix = getMatrix(x, y, 0)
	traceBack = TraceBack(x, y)

	for i in range(1, len(x) + 1):
		for j in range(1, len(y) + 1):
			hor = matrix[i][j - 1] + score.gap
			if hor < 0:
				hor = 0
			ver = matrix[i - 1][j] + score.gap
			if ver < 0:
				ver = 0
			diag = matrix[i - 1][j - 1] + score.misMatch(x[i - 1], y[j - 1])
			if diag < 0:
				diag = 0
			matrix[i][j] = max(hor, ver, diag)
			if matrix[i][j] == hor:
				traceBack[i][j] = 'hor'
			elif matrix[i][j] == ver:
				traceBack[i][j] = 'ver'
			else:
				traceBack[i][j] = 'diag'
	return matrix, traceBack


###################### MAIN CODE ##########################

fileName1 = input("Enter first sequence file: ")
file1 = open("./P1AASeqs/" + fileName1)
sequence_1 = file1.readlines()[1]

fileName2 = input('Enter second sequence file: ')
file2 = open("./P1AASeqs/" + fileName2)
sequence_2 = file2.readlines()[1]

fileNameSUB = input("Enter sub matrix file: ")
input1 = np.loadtxt("./P1SubMatrices/" + fileNameSUB, skiprows=2, dtype='i', delimiter=',')
print(input1)

gapPenalties = input("\nEnter Gap Penalties: ")

print("\nMatch and MisMatch values are supposed to set by the sub matrix")
match = input("\nEnter Match Value: ")
misMatch = input("Enter MisMatch value: ")

x = sequence_1
y = sequence_2

print('\nInput sequences are: ')
print(x)
print(y)
print()

score = Scoring(int(gapPenalties),int(match),int(misMatch))

mode = input("Global, Semi-Global, Local :")
if mode == "Global":
	matrix,traceBack = globalAlign(x,y,score)
elif mode == "Semi-Global":
	matrix, traceBack = semiGlobalAlign(x, y, score)
else:
	matrix, traceBack = localAlign(x, y, score)

print('\nScore matrix:')
for i in range(len(matrix)):
	print(matrix[i])
print()

xSeq,ySeq = getAlignedSequences(x,y,matrix,traceBack)

print('\nThe globally aligned sequences are:')
print(*xSeq[::-1])
print(*ySeq[::-1])