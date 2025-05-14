############################################################################
# Connor Bingham
# CutsiteFinder
# Finds a given RNA sequence in a single or multiple genes (in FASTA format)
# If given a secondary structure beneath RNA sequence, will also reveal if the sequence lies on a hairpin loop
# Created: 6/26/2024
# Last Updated: 5/6/2025
############################################################################
import sys 

def CreateSecStruct(dotBracket,SECONDARY_STRUCTURES): #translates the dot-bracket list into what the secondary structure actually means. This does not support Pseudoknots or G-Coupling
    finalSecStruct = []    
    curStruct = SECONDARY_STRUCTURES[0]
    oldStruct = SECONDARY_STRUCTURES[2]
    loopLayers = 0 #this might be a useful variable if you are seeking to expand the program to search for more complex structures; it increases each time the structure "branches" and decreases each time it "debranches."
    lastStem = 0
    isOnLoop = False
    isThreePrimeStem = False
        
    for nucleotide in range(len(dotBracket)):
        if dotBracket[nucleotide] == "(" or dotBracket[nucleotide] == ")": #handling stems
            curStruct = SECONDARY_STRUCTURES[1]
            
            if dotBracket[nucleotide] == "(":
                loopLayers += 1
                isThreePrimeStem = False
            else:
                loopLayers += -1
                isThreePrimeStem = True
            
        elif dotBracket[nucleotide] == ".": #handling unpaired regions, and declares if a nucleotide is on any form of loop
            if loopLayers <= 0:
                curStruct = SECONDARY_STRUCTURES[0]
            else:
                if isThreePrimeStem:
                    curStruct = SECONDARY_STRUCTURES[2]
                else:
                    if curStruct == SECONDARY_STRUCTURES[1]:
                        lastStem = nucleotide
                    isOnLoop = True
                    curStruct = "None"
        
        if curStruct == SECONDARY_STRUCTURES[1] and isOnLoop == True: #determines if a previous loop was a hairpin loop or a interior loop. 
            if isThreePrimeStem:
                oldStruct = SECONDARY_STRUCTURES[3]
            else:
                oldStruct = SECONDARY_STRUCTURES[2]
            count = lastStem
            while count < nucleotide:
                finalSecStruct[count] = oldStruct
                count += 1
            isOnLoop = False
        
        finalSecStruct.append(curStruct)
    return finalSecStruct

def Findgene(sequence,file,name): #performs the actual data anaylsis, given a gene with a dot-bracket secondary structure
    data = Buildgene(file)
    gene = data[0].upper()
    
    ##### This translates the dotbracket into its actual meaning
    SECONDARY_STRUCTURES = ("Unpaired","Stem","Interior Loop","Hairpin Loop") #this program is not build to find multiloops/junctions. Unpaired refers to incomplete loops that can be at the very start and end of mRNAs
    secStruct = CreateSecStruct(data[1],SECONDARY_STRUCTURES)
    assert len(secStruct) == len(gene), input(name+" gene's secondary structure does not match the sequence length, please correct this and try again")
    #####
    
    startline = 0
    position = 0
    onSeq = 0
    matchLoc = []
    genLen = len(gene)
    nucleotide = 0 #treat this identically to a count
        
    hairpinStart = 0
    curHairpinSeq = ""
        
    while nucleotide < genLen:
        
        #The program is meant to locate only hairpin-loop structures, but can be modified to search for other structures by changing every instance of SECONDARY_STRUCTURES[3] to a different value in the list. This hasn't been tested effectively, so there may be unexpected bugs. The printing portion of the code will need to be updated seperately if that is something you deem important.
        if secStruct[nucleotide] != SECONDARY_STRUCTURES[3]: 
            curHairpinSeq = ""
        else:
            #This gets the extent of what the complete area of the local structure is. Important to note, if you switch this to not using hairpin loops anymore, this will not be accurate, and will be missing later and earlier nucleotides (for example, if you switched it to searching for stem structures, it would only return whatever half of the stem the sequence is on, missing the other half)
            if curHairpinSeq == "":
                hairpinStart = nucleotide
                while secStruct[nucleotide] == SECONDARY_STRUCTURES[3]: 
                    curHairpinSeq += gene[nucleotide]
                    nucleotide += 1
                nucleotide = hairpinStart
            
        ######################################
        
        #This section matches up the sequence you asked it to find, to the gene sequence of an mRNA, to see if there are any matches. The "elif" checks to see if the your first letter matches the current nucleotide the program is looking at, and the "if" checks whether it continues to align with the sequence.      
        if onSeq > 0:
            if gene[nucleotide] in sequence[onSeq]:
                onSeq += 1
                if not secStruct[nucleotide] in seqSecStruct:
                    seqSecStruct.append(secStruct[nucleotide])
                        
                if onSeq >= len(sequence):
                    specificSeq = gene[seqStartLoc:nucleotide+1]  #this gives what specific sequence . E.x., if sequence was GNG, it finds what nucleotide N is in this case
                        
                    listForEachInstanceOfSeqeunce = []
                    listForEachInstanceOfSeqeunce.append(seqStartLoc+2)
                    localSeq = gene[seqStartLoc-4:seqStartLoc].lower()+specificSeq+gene[nucleotide+1:nucleotide+5].lower()
                    listForEachInstanceOfSeqeunce.append(localSeq)
                    listForEachInstanceOfSeqeunce.append(seqSecStruct)
                        
                    ############################ 
                    if seqSecStruct == [SECONDARY_STRUCTURES[3]]:
                        localHairpin = curHairpinSeq[0:seqStartLoc-hairpinStart].lower()+specificSeq+curHairpinSeq[nucleotide+1-hairpinStart:len(curHairpinSeq)].lower()
                        listForEachInstanceOfSeqeunce.append(localHairpin)
                    ############################
                            
                    matchLoc.append(listForEachInstanceOfSeqeunce)
                        
                    nucleotide = nucleotide - (onSeq-1)
                    onSeq = 0
            else:
                nucleotide = nucleotide - (onSeq-1)
                onSeq = 0
            
        elif gene[nucleotide] in sequence[0]:
            seqStartLoc = nucleotide
            seqSecStruct = [secStruct[nucleotide]]
            onSeq = 1
        nucleotide += 1
    return matchLoc
    
    
def Buildgene(gene): #seperates a FASTA file into a gene sequence and a correspoding secondary structure (it is has both). This is built with the assumption that the sequence comes first, the secondary structure comes second. 
    baseGene = gene
    baseGene.readline()
    
    finalGene = ""
    finalDotBracket = ""
    loop = True
    contents = "gene"
    while loop:
        nextLine = baseGene.readline()
        if "(" in nextLine or ")" in nextLine or "." in nextLine or "1" in nextLine: #this determines whether the section is a dot-bracket or not. 
            contents = "DotBracket"
            
        if contents == "gene":
            if "U" in nextLine.upper(): # Convert mRNA to cDNA 
                nucleotide = 0
                nextLine = list(nextLine)
                for char in nextLine:
                    
                    if char.upper() == "U":
                        nextLine[nucleotide] = "T"
                    nucleotide += 1
                nextLine = "".join(nextLine)
            finalGene += nextLine.strip("\n")
        elif contents == "DotBracket":
            finalDotBracket += nextLine.strip("0123456789 \n")
            
        if nextLine == "":
            loop = False

    return (finalGene, finalDotBracket)
    
  
def genHub(wholeList,sequence): #is responsible for finding each of the files that are in the Genlist, and then running the anaylsis functions on them
    matchList = []
    for gene in wholeList:
        seqList = []
        gene = gene.strip("\n")
        seqList.append(gene) 
        
        try:
                genFile = open("Genes\\"+gene+".txt","r")
        except: 
                input(gene+".txt"+" is not a file in the gene folder. Make sure it is capitalized properly and in the gene folder.")
                print(intentionalCrash)          
        
        seqList.append(Findgene(sequence,genFile,gene))
        genFile.close()
        
        print(gene,"complete",file=sys.stderr)
        matchList.append(seqList)
    return matchList
    

def writeResults(matchList,sequence,showOnlyYes): #writes the results of the experiment to stdout
    matches = 0
    matchesOver1 = 0
    hairPins = 0
    partialHairPins = 0
    lines = ""
    for gene in matchList:
        genName = gene[0]
        
        if not gene[1]:
            lines += genName+" --> No matches"+"\n" 
        else:
            matches += 1
            if len(gene[1]) > 3:
                matchesOver1 += 1
            
            lines += genName+" --> "+str(len(gene[1]))+" matches"+"\n"
            
            count = 1
            genHairPins = 0
            genPartialHairPins = 0
            
            
            for match in gene[1]:
                tempLine = "Match "+str(count)+" is at position "+str(match[0]-1)+"\n"
                count += 1
                    
                if len(match) > 1:
                    tempLine += "sequence:  "+match[1]+"\n"
                
                if len(match) > 2:
                    if match[2] == ["Hairpin Loop"]: 
                        isHairpin = "Yes"
                        genHairPins = 1
                        genPartialHairPins = 1
                        
                    elif "Hairpin Loop" in match[2]:
                        isHairpin = "Partially"
                        genPartialHairPins = 1
                        
                    else:
                        isHairpin = "No"
                        
                    tempLine += "is on hairpin-loop?:  "+isHairpin+"\n"
                if len(match) > 3:
                    tempLine += "Sequence of hairpin-loop: "+match[3]+"\n"
                    tempLine += str((len(match[3])-len(sequence)))+" additional nucleotides on the loop \n"
                
                if not showOnlyYes or isHairpin == "Yes":                
                    lines += tempLine
                
            hairPins += genHairPins
            partialHairPins += genPartialHairPins
        lines += "\n"
    print("Searched for sequence:",sequence)
    print(len(matchList),"Genes were searched")
    print()
    print(matches,"Genes have the Sequence")
    print(len(matchList)-matches,"Genes do not have the Sequence")
    print(matchesOver1,"Genes have more than one instance of the Sequence")
    print()
    print(hairPins,"Genes have the sequence entirely on a hairpin-loop")
    print(partialHairPins,"Genes have the sequence at least partially on a hairpin-loop")
    print()
    percent = matches/len(matchList)*100
    print("%0.2f Percent of genes have a %s sequence" %(percent,sequence))
    percent = matchesOver1/len(matchList)*100
    print("%0.2f Percent of genes have more than one %s sequence" %(percent,sequence))
    percent = hairPins/len(matchList)*100
    print("%0.2f Percent of genes have a %s sequence on a hairpin-loop" %(percent,sequence))
    percent = partialHairPins/len(matchList)*100
    print("%0.2f Percent of genes have a %s sequence at least partially on a hairpin-loop" %(percent,sequence))    
    
    print()
    print("-------------------------")
    print(lines)
  
    
def generateNucleotides(): #this allows nucleotides like K M and N to be properly recognized and dicypher what nucleotides they represent. 
    nucFile = open("Parameters\\"+"Valid Nucleotides.txt","r")
    nucListRaw = nucFile.readlines()
    count = 0
    nucList = []
    
    for nuc in nucListRaw:
        nuc = nuc.strip("\n")
        nuc = nuc.split()
        nucList.append(nuc)
    
    nucFile.close()
    return nucList


def main():
    print("Please input the name of the list containing the names of genes to search", file=sys.stderr)
    listName = input()
    
    if not ".txt" in listName:
        listName += ".txt"
        
    try:
            readingList = open("Parameters\\"+listName,"r")
    except: 
            assert 1==2,input("The file you have requested does not exist, please restart the program and try again")
    
    sequence = False
    nucReference = generateNucleotides()
    
    validNuc = []
    for nuc in nucReference: #this generates a list of valid letters for inputed sequences
        validNuc.append(nuc[0])

    while not sequence:
        print("Please input the sequence to be found", file=sys.stderr)
        sequence = input().upper()
        isValid = True
        rebuiltSeq = ""
        
        if len(sequence) <= 1:
            isValid = False
            print("\nSequences must be at least two nucleotides long", file=sys.stderr)          
        
        for letter in range(len(sequence)):
            if sequence[letter] not in validNuc:
                isValid = False
                print("\nSequence invalid, please enter a valid series of nucleotides", file=sys.stderr)                 
            elif sequence[letter] == "U":
                rebuiltSeq += "T"
            else:
                rebuiltSeq += sequence[letter]
                      
        if not isValid:
            sequence = False
        else:
            sequence = rebuiltSeq
            
    finSequence = []
    for letter in sequence: #this converts the sequence to its corresponding AGCT nucleotides
        for nuc in nucReference:
            if letter == nuc[0]:
                finSequence.append(nuc[1])
            
    assert len(finSequence) == len(sequence), input("An error occured when converting inputed sequence to corresponding nucleotides, please report this bug")
    
    print("Only display results with matching secondary structures?", file=sys.stderr)
    showOnlyYes = input().upper()    
    showOnlyYes = showOnlyYes == "Y" or showOnlyYes == "YES"
        
    print("", file=sys.stderr)
    print("", file=sys.stderr)
    writeResults(genHub(readingList.readlines(),finSequence),sequence,showOnlyYes)
    
    readingList.close()
    print("", file=sys.stderr)
    print("Input anything to end the program", file=sys.stderr)
    input()
    
main()