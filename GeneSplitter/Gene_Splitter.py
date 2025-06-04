############################################################################
# Made by: Connor Bingham
# Gene_Splitter
# Made to turn a single text file with a large list of genes into individual files for each gene. Also works with output (as a text file) from ViennaRNA 
# Created on: 7/19/2024
# Last Updated: 5/6/2025
############################################################################
import sys 

print("Please input the name of the master file", file=sys.stderr) #obtain files and suffix from users
listName = input()
if not "." in listName:
    listName += ".txt"
try:
    genFile = open(listName,"r")
except: 
    print("The file you have requested does not exist, please restart the program and try again", file=sys.stderr)
    input()
    intentionalCrash
        
print("Enter a suffix to file names (to indicate species); leave blank if you want no suffix", file=sys.stderr)
suffix = input()

loop = True
inFile = False
genList = ""

while loop:
        nextLine = genFile.readline()
        if nextLine == "": #close once the list is over
                loop = False   
                if inFile:
                        gene.close()
        elif ">" == nextLine[0]: #Handle New Gene
                if inFile:
                        gene.close()
                        
                if listName[len(listName)-4:len(listName)-4] == ".fna":
                    genName = nextLine.split()
                    genName = genName[1]
                else:
                    genName = nextLine.strip(", mRNA\n").strip(">")
                    
                    testForFASTA = genName.split()
                    if "(" in testForFASTA[len(testForFASTA)-1] and ")" in testForFASTA[len(testForFASTA)-1]:
                        genName = genName.split()
                        genName = genName[len(genName)-1].strip("(").strip(")")
                    
                print(genName+suffix, file=sys.stderr)
                genList += genName+suffix+"\n"
                
                gene = open("SplitList/"+genName+suffix+".txt","w")
                gene.write(nextLine)
                inFile = True
        else:
                if inFile: #write Gene Contents 
                    if "(" in nextLine or "." in nextLine or ")" in nextLine: #checks if this is a secondary structure or a 
                        cleanNextLine = nextLine.split()
                        nextLine = cleanNextLine[0]
                    
                    gene.write(nextLine)
    
genListFile = open("Genlist"+suffix+".txt","w") #create a Genlist for CutsiteFinder.py
genListFile.write(genList) 
genListFile.close()
                    
genFile.close()