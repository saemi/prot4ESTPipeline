##########################################################################################################
###This script reads in protein sequence from prot4Test and 
###1)removes short proteins (default shorter than 100 aa)
###2)Uses the cluster output and retrives the 'best' cluster 
###3)removes stop codons and undefined amino acids (i.e. X)
###4)Renames the output contigs for orthoMLC (keeps a record of what contigs were used)
###5)Retrives the nucl format 
##########################################################################################################
###
###THINGS TO ADD & FIX
### 1. MAKE NUCLEOTIDE SEQUENCE COME IN SAME ORDER AS PROTEIN SEQUENCES
### 2. MAKE THE STOP CODON PENALTY A VARIABLE, WITH A DEFAULT OF NEGATIVE 10
##########################################################################################################
from Bio import SeqIO
import sys, fnmatch

inputfastaProt=sys.argv[1]
#inputfastaProt = 'L.bienne_p4ePro.fsa'

inputfastaNucl = sys.argv[2]
#inputfastaNucl = 'L.bienne_p4eNuc.fsa'
###Uncomment if you want to have it automatic
#inputfastaNucl = str(inputfastaProt).split('_')[0]+'_p4eNuc.fsa'

###Here I could make the script pull it from the from the contig names, see line ~180
speciesAbbrevation = str(sys.argv[3])
#speciesAbbrevation = 'Lbie'

###The sequence length is optional, but if you don't put a value, it puts 100aa
try: 
	sequenceLength = int(sys.argv[4])
except IndexError:
	sequenceLength = 100
	pass

##########################################################################################################

def removeSmallTranscripts(sequences):
	'''This function takes sequences SeqRecord format and returns a dictonary
	with seq ID as keys and sequences as values'''
	largeContigName = []
	largeContigSequence = []
	sequences = recordsProt
	for i in recordsProt:
		if len(i) > sequenceLength-1:
			largeContigName.append(removeUnderScore(i.id))
			largeContigSequence.append(i.seq)
	contigDict = dict(zip(largeContigName,largeContigSequence))
	return contigDict	


def removeUnderScore(contigID):
	'''prot4EST gives contigs identified with blast an _1 so I just remove
	it here but it might be worthwhile at some point to do something with that info'''
	return contigID.replace("_1","")


def groupClusteredTranscripts(contigNames):
	'''This funcitons identifes transcripts that are in clusters and puts them in a dictionary'''
	largeContigClusterName = []
	largeContigClusterSequence =[]
	contigNames = contigDict.keys()
	#for i in largeContigName:
	for i in contigNames:
		if int(i.split('.')[1]) > 1:
			try:
				contigDict[str(i.split('.')[0] + '.' + str(1))]
				largeContigClusterName.append(str(i.split('.')[0] + '.' + str(1)))
				largeContigClusterSequence.append(contigDict[str(i.split('.')[0] + '.' + str(1))])
			except KeyError:
				pass		
			for h in range(2, int(int(max(sizeOfCluster))+1)):
				try:
					contigDict[str(i.split('.')[0] + '.' + str(h))]
					largeContigClusterName.append( str(i.split('.')[0] + '.' + str(h)))
					largeContigClusterSequence.append(contigDict[str(i.split('.')[0] + '.' + str(h))])
				except KeyError:
					pass
	clusterContigDict = dict(zip(largeContigClusterName,largeContigClusterSequence))
	return clusterContigDict


def sortClusteredTranscripts(clusterContigDict):
	'''Sorts the clusters and returns a tuple of two lists. The output is used in
	the translateClusteredTranscripts function'''
	largeContigClusterNameSorted = []
	largeContigClusterSequenceSorted = []
	for key in sorted(clusterContigDict.iterkeys()):
		largeContigClusterNameSorted.append(key)
		largeContigClusterSequenceSorted.append(clusterContigDict[key])
	return zip(largeContigClusterNameSorted,largeContigClusterSequenceSorted)


def removeXandStar(proteinSequence):
	'''when there is only one sequence per cluster, this removes stars and X aa'''
	seq = str(proteinSequence)
	#return seq.replace("X","").replace("*","")
	###OrthoMLC might want to have the *, to estimate protein quality
	return seq.replace("X","")
	

###Here might be good to report the name of full name, with the dot(e.g. .3, .4 or .5)
###Do that for the next version of the script
def translateClusteredTranscripts(sortedClusterOfContigs):
	'''This function groups transcript clusters ands sends them to be tested
	for transcript quality.'''
	x=0
	bestClusteredContigName = []
	bestClusteredContigSequence = []
	largeContigClusterNameSorted = []
	largeContigClusterSequenceSorted = []
	for i in sortedClusterContigs:
		largeContigClusterNameSorted.append(i[0])
		largeContigClusterSequenceSorted.append(i[1])
	for contigs in largeContigClusterNameSorted:
		while len(largeContigClusterNameSorted) > 0:
			contigOfInterest=largeContigClusterNameSorted[x].split('.')[0]+'.'+ '*'
			listOfClusteredContigs = fnmatch.filter(largeContigClusterNameSorted,str(contigOfInterest))
			#print '>' + contigOfInterest.split('.*')[0] + '\n' + str(translateBestTranscript(listOfClusteredContigs))
			bestClusteredContigName.append(contigOfInterest.split('.*')[0])
			bestClusteredContigSequence.append(str(translateBestTranscript(listOfClusteredContigs)))
			for i in listOfClusteredContigs:
				largeContigClusterNameSorted.pop(largeContigClusterNameSorted.index(i))
	return dict(zip(bestClusteredContigName,bestClusteredContigSequence))


###CHANGE THE VALUES HERE YOU WANT TO ADJUST THE CRITERA OF TRANSCRIPT SELECTION
###Also I could make my internal stop codon penalty a varialbe 
def transcriptQuality(nameOfContig):
	'''This function calculates the qulaity of the transcript and takes
	the length of the transcript, plus the number and position of stop codons into account'''
	proteinSeq = str(clusterContigDict[nameOfContig])
	###When no stop codons
	if proteinSeq.count('*') == 0:
		quality = len(proteinSeq)	
	else:
		###If only one at the end, does not matter
		if proteinSeq.find('*') == 1-len(proteinSeq):
			quality = len(proteinSeq)
		elif proteinSeq.rfind('*') == 1-len(proteinSeq):
			quality = (len(proteinSeq) - 10*proteinSeq.count('*'))-10
		else:
		###But all internal stop codons give a score of -10 
			quality = len(proteinSeq) - 10*proteinSeq.count('*')
	return quality


def translateBestTranscript(ClusteredContigs):
	'''This identifies the highest transcripton quality and returns a dictonary
	with contig names and cleaned protein sequences'''
	bestClusteredContigs = []
	bestClusteredContigsName = []
	bestClusteredContigsSequence = []
	contigsInQuestionName = []
	contigsInQuestionSequence = []
	contigsInQuestionQuality = []	
	if len(ClusteredContigs) == 1:
		bestClusteredContigs.append(ClusteredContigs)
		return bestClusteredContigsSequence.append(clusterContigDict[ClusteredContigs[0]])
	else:
		lengthOfCluster = len(ClusteredContigs)
		for contig in ClusteredContigs:
			contigsInQuestionQuality = []
			contigsInQuestionName.append(contig)
			contigsInQuestionSequence.append(clusterContigDict[contig])
			for contigInQuestion in contigsInQuestionName:
				contigsInQuestionQuality.append(transcriptQuality(contigInQuestion))
		bestContigIndex = contigsInQuestionQuality.index(max(contigsInQuestionQuality))
		return clusterContigDict[contigsInQuestionName[bestContigIndex]]


###Need to fix this so it keeps the orignal name with all of the information, not just id
###Well I could do that for the next, I don't really need it 
def contigNameGenerator(contigName):
	'''Renames the contigs for otrthoMCL and writes out a list of names used'''
	oldName = str(contigName)
	oldNameList.append(oldName)
	newName =  '>' + speciesAbbrevation + '|' + str(numberOfContigs.pop(0))
	newNameList.append(newName)
	return newName


def findTheBiggestCluster(recordsProt):
	'''Makes a list of the sizes of clusters'''
	sizeOfCluster = []
	for contig in recordsProt:
		sizeOfCluster.append((removeUnderScore(contig.id)).split('.')[1])
	return sizeOfCluster

def removeIdenticalProteins(contigDictNoClusters,contigDictWithBestClusters):
	'''This fucntion removes identical protins from the dictionaries, by combining two dictonaries
	and swapping the values and keys once, thereby removing non-uniqe values by using them as keys'''
	combinedDict = dict(contigDictNoClusters.items() + contigDictWithBestClusters.items())
	combinedDict = dict(contigDict.items() + bestClusteredContigDict.items())
	combinedDictNoDupsName = []
	combinedDictNoDupsSequence = []
	for contig in combinedDict.keys():
		combinedDictNoDupsName.append(contig)
		combinedDictNoDupsSequence.append(removeXandStar(combinedDict[contig]))
	combinedDictNoDups=dict(zip(combinedDictNoDupsSequence,combinedDictNoDupsName))
	#combinedDictNoDupsReturn=dict(zip(combinedDictNoDups.values(), combinedDictNoDups.keys()))
	return dict(zip(combinedDictNoDups.values(), combinedDictNoDups.keys()))

def parseFastaNucl(recordsNucl):
	'''This function gets the nuclotide sequences of the protein translated
	transcript that are pulled out of inputfastaProt'''
	tmpNuclContigName = []
	tmpNuclContigSequence = []
	nuclContigName = []
	nuclContigSequence = []
	for i in recordsNucl:
		tmpNuclContigName.append(removeUnderScore(i.id))
		tmpNuclContigSequence.append(i.seq)
	tmpDict = dict(zip(tmpNuclContigName,tmpNuclContigSequence))
	for nuclContig in recordsNucl:
		if removeUnderScore(nuclContig.id) in combinedNonDuplicatedProteinDict.keys():
			nuclContigName.append(removeUnderScore(nuclContig.id))
			nuclContigSequence.append(tmpDict[removeUnderScore(nuclContig.id)])
		elif removeUnderScore(nuclContig.id).split('.')[0] in combinedNonDuplicatedProteinDict.keys():
			nuclContigName.append(removeUnderScore(nuclContig.id).split('.')[0])
			nuclContigSequence.append(tmpDict[removeUnderScore(nuclContig.id)])
	return dict(zip(nuclContigName,nuclContigSequence))
	
##########################################################################################################

handle = open(inputfastaProt, "rU")
recordsProt = list(SeqIO.parse(handle, "fasta"))
handle.close()

handle = open(inputfastaNucl, "rU")
recordsNucl = list(SeqIO.parse(handle, "fasta"))
handle.close()

#speciesAbbrevation = (records[0].id).split('_')[0]
sizeOfCluster = findTheBiggestCluster(recordsProt)

###These dictonaries keep different versions of the transcripts
contigDict = removeSmallTranscripts(recordsProt)
clusterContigDict = groupClusteredTranscripts(contigDict.keys())
sortedClusterContigs = sortClusteredTranscripts(clusterContigDict)
bestClusteredContigDict = translateClusteredTranscripts(clusterContigDict)
combinedNonDuplicatedProteinDict = removeIdenticalProteins(clusterContigDict,bestClusteredContigDict)

###Used to name contigs and keep them in two lists
numberOfContigs = range(1,(len(bestClusteredContigDict)+len(contigDict)+1))
oldNameList = []
newNameList = []

###remove clustered contigs from contig dictionary
for i in clusterContigDict.keys():
	contigDict.pop(i)
	
###The Outputs
###If I want information on what names where used, add the write out function
outputfileProt = speciesAbbrevation + '.UnigenesProt4EST_prot.fsa'
outputfileNucl = speciesAbbrevation + '.UnigenesProt4EST_nucl.fsa'
outputfileContigNames = speciesAbbrevation + '.UnigenesProt4EST_namesUsed'

###This gets the contigs in nucleotide format
nuclContigDict = parseFastaNucl(recordsNucl)

numberOfContigs = range(1,(len(combinedNonDuplicatedProteinDict)+1))
f1 = open (outputfileProt,"w")
for i in combinedNonDuplicatedProteinDict.keys():
	f1.write(contigNameGenerator(i) + '\n' + str(combinedNonDuplicatedProteinDict[i])+'\n')
f1.close()

###This dictonary has to be here because the contigNameGenerator function makes the names
newNameOldNameDict = dict(zip(oldNameList,newNameList))

f2 = open (outputfileNucl,"w")
for i in nuclContigDict.keys():
	f2.write(newNameOldNameDict[i] + '\n' + str(nuclContigDict[i])+'\n')
	#print newNameOldNameDict[i], nuclContigDict[i]

f2.close()

f3 = open (outputfileContigNames,"w")
for i in zip(newNameList,oldNameList):
	f3.write(i[0]+'\t'+'>'+i[1]+'\n')
f3.close()

##########################################################################################################
###Uncomment these if you don't want to exclue non unique proteins
#f1 = open (outputfileProt,"w")
#for i in contigDict.keys():
	#f1.write(contigNameGenerator(i) + '\n' + removeXandStar(str(contigDict[i]))+'\n')
	#
#f1.close()
#
#
#f1 = open (outputfileProt,"a")
#for i in bestClusteredContigDict.keys():
	#f1.write(contigNameGenerator(i) + '\n' + removeXandStar(str(bestClusteredContigDict[i]))+ '\n')
#
#f1.close()
