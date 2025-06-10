# CutSequenceFinder

A package that can search specified genes for specified sequences located on a specified structure. 

----

## Table of Contents
1. [Prerequisite Knowledge](#Prerequisites)
2. [The SeqSearch Program](#SeqSearch)
3. [Obtaining Secondary Structures](#Secondary)
4. [Reference Guide](#Reference)
5. [License](#License)
6. [Contact](#Contact)

----

## Prerequisites

This package assumes that you have installed Python 3.0 and know how to run Python programs from the terminal, and how to do basic input and output redirection. If you do not, there are plenty of resources online to help you with these components. (note that this is different on Windows and Linux systems)
It is also assumed that you have gene sequences that you plan to search within or know how to obtain them. (If this is not true, [Genbank][Genbank] usually hosts most genes) 

## SeqSearch

The SeqSearch Program is the primary program of the CutSequenceFinder, and can be used independently from the GenSplitter program. SeqSearch is responsible for actually doing the searching through genes with secondary structures for a specific sequence.

1. Structure of SeqSearch
    
    SeqSearch contains three important components: the Genes Folder, the Parameters folder, and the SeqSearch Python script. 

    The Genes folder should contain all of the genes that can be searched (either as RNA or cDNA). When you have genes that you want to search, make sure they are in this folder. All genes in the folder should have a dot-bracket secondary structure included. How to obtain this will be mentioned [later](#Secondary Structures)

    The Parameters folder is how you specify what specific genes in the Genes Folder will be searched for any individual run of SeqSearch. To facilitate this, you can make a text file with the name of each gene that you want to search for (one gene per line) in this folder. The name in the text file must match the name in the genes folder perfectly. 
    If you include genes in the text file that are not in the Genes Folder, the program will not work. 
    There is also a file in the Parameters folder named "Valid Nucleotides.txt", which cannot be removed. This folder enables support for searching using all IUPAC nucleotide codes. 

    Finally, there is the SeqSearch Python script, which can be executed to run the program. There are no executables associated with this program. 

2. Using SeqSearch
    
    To use the SeqSearch program, run the SeqSearch python script. The python script itself contains some instructions for using the program, that will be repeated and expanded on here. 

    When run, the SeqSearch program needs specification of:
    - A Text File in the parameters folder (specifying what genes will be searched)
    - What sequence should be searched for
    - Whether you want to see all instances of the sequence or just the ones on the proper structure (yes or no)

    The details of the text file are mentioned in step 1. If you do not specify the extension of the file (such as Genlist instead of Genlist.txt), it will automatically assume that you are refering to a .txt file. 

    For the sequence to be searched it must be at least 2 nucleotides long, and only contain valid nucleotides. There is no functional difference between T and U for this program, and the extended nucleotide dictionary (such as N for any nucleotide, and R for the purines) is supported. 

    It is worth mentioning that the yes or no specification is actually implemented by checking whether the first letter is Y or not. That means inputting anything that starts with Y will be interpreted as yes, and anything that does not start with Y will be interpreted as no. 
    What this actually does is not visually declare instances of the sequence that aren't on the proper secondary structure. This is purely to reduce clutter if you only care about matches that are on the desired secondary structure. More details are in step 3. 

    To record the output; for runs of only a few genes (<100), you can simply copy and paste the terminal output into somewhere more permanent. For larger runs (>100), it is recommended that you redirect the output of the SeqSearch program to a text file.
    (The unimportant output (The instructions on how to use the program and progress tracking) is printed to standard error instead of standard output, so they will automatically not be included if you redirect the output of SeqSearch) 

3. Output
    
    There are two components to the output of this program; the summary and the individual gene results. They are separated by a series of dashes in the output.  

    The summary contains a helpful depiction of the results across all of the genes searched. Note that "genes having the sequence" refers to the gene having any instance of the sequence at all, not that the sequence is specifically on the proper secondary structure. The results for that are specified directly.

    For the individual gene results, the results of each gene are separated by a paragraph. The individual gene results are less obvious, so the structure is as follows:

    The first line contains the name of the gene and how many instances of the sequence were contained in it. Again, this refers to any match, not matches on a specific secondary structure. 
    Then each match will be listed (in order of first appearance), with details on the match specified. If you specified "yes" to only displaying results with proper secondary structures, any matches that are not on the proper secondary structure will not be displayed. 

    The first line of each match describes the position of the sequence on the gene (in nucleotides).
    The second line of each match shows the sequence and the surrounding nucleotides (to make it easier to find in the gene)
    The third line of each match specifies whether the match is on the proper secondary structure (If it is on multiple structures including the proper structure, it will specify 'partially') The next lines are only present if the answer is yes. 
    
    The fourth line of each match will show all nucleotides on the loop that the sequence is located on (this only works with hairpin loops)
    The fifth line of each match will say the length of the loop the sequence is located on (this only works with hairpin loops)

4. Limitations

    SeqSearch by default can only target structures on hairpin loops, and cannot directly search for sequences on other RNA secondary structures. The SeqSearch program can be fairly easily modified to search for other secondary (instructions for which are contained in the python script itself); but some parts of the output will not be returned properly (such as incorrectly naming the output as "on Hairpin-loops"). Found Matches and locations of the found Matches will always be correct (unless the program is modified incorrectly)

    SeqSearch also only supports the basic notation of the Dot-Bracket Notation. (only dots and parentheses) There is no support for Psuedoknots or G-coupling, for example.

    The program would also require more extensive modification to distinguish junctions/multiloops from normal internal loops, so it is not recommended to use this program for that purpose. 

## Secondary Structures

To use SeqSearch, all genes in the gene folder need to have a name in the first line, then a sequence, and then the secondary structure for that sequence in dot-bracket notation. This section details a few ways to get genes with secondary structures into the program.

1. ViennaRNA

    ViennaRNA is not the only way to obtain secondary structures, and different ways should work as long as you obtain a properly formatted gene. (Although only ViennaRNA has been tested) Details for reference are described on the ViennaRNA [website][vrna_website].

    (One-by-one)
    If doing genes one at a time, or using a small number of genes, it is recommended to use the [RNAfold][RNAfold] web server hosted on the ViennaRNA website for simplicity. How to use RNAfold is described on the website

    Relevent to this program is that the website will return a secondary structure in the form of a dot-bracket notation. This is the part with the parentheses and periods below your gene sequence.
    (It will have two different secondary structures; which refer to different prediction methods. Either is a valid secondary structure, but one may be more accurate than the other)

    For later purposes, you should make a text file with the name of your gene as the name of the file and put the sequence of the gene into the file first, and then copy and paste the dot-bracket notation from ViennaRNA into the program second. 
    (It is okay to copy the numbers listing the sequence position when you copy the secondary structure, the program automatically removes them)

    (in bulk)
    If doing in bulk, first make sure to obtain the [ViennaRNA Package][vrna_package] for your operating software and understand how to access and use the RNAfold program. You do not have to compile the binary package from source code, and the executables for Windows and MacOS X will work.  

    To obtain a file containing proper secondary structures, use RNAfold and redirect the output of the process to a text file. This text file does not need to be named anything in particular. Make sure the genes that you put into ViennaRNA are in FASTA format; where the first line is the name of the gene (indicated by >), and the sequence comes after. This is not necessary for ViennaRNA to work, but it is helpful later in the process. 

    If you are working with a large file containing multiple genes (such as Genome files from Genbank), it is recommended to add the "@" character to the end of the file and redirect the input of RNAfold to that file. (the @ character is added so that RNAfold stops once it is complete).

    Note: If the name of your gene contains certain special characters (such as /), RNAfold will not stop processing genes once it arrives at that character. Make sure to change these characters if you encounter problems. 


2. Importing to CutSequenceFinder 

    As mentioned previously, for genes to be searched by CutSequenceFinder, they must be placed in the Genes folder of SeqSearch. This section details how that should be done. (You need to already have secondary structures at this point)
    
    (One gene per file)
    For importing genes where each gene is already its own file, you can directly move the text file previously mentioned into the genes folder in SeqSearch. Although, do make sure it is a .txt text file. It will not work otherwise. 

    (Multiple genes in a file)
    For importing genes that have multiple genes in one file (such as the redirected output of RNAfold), CutSequenceFinder provides a program to aid in transferring these genes to the genes folder, GeneSplitter. GeneSplitter is only designed to work with files in FASTA format (with each gene's secondary structure being after its sequence), so ensure that your file is in that format.  

    To use the GeneSplitter program, place your file into the GeneSplitter folder, and run the GeneSplitter Python script. GeneSplitter only needs you to specify the name of your file to function (adding a suffix is optional). If you do not specify the extension of the file (such as Genome instead of Genome.txt), it will automatically assume that you are referring to a .txt file.

    When run, GeneSplitter will make a file for each individual gene and place them in the genes folder of SeqSearch. The name of each file will be its FASTA header, plus the suffix (if you specify one). GeneSplitter does not check to see if the file already exists in the genes folder, so files with the same name will be overwritten. 
    Like RNAfold, if the name of a gene contains certain special characters (such as /), GeneSplitter will not be able to process it. Make sure to change these characters if you encounter problems.

    GeneSplitter will also automatically create a text file in the parameters folder containing the name of every gene it found, that can be used in the SeqSearch program. This is by default called "Genlist(Suffix).txt", so it is recommended that you change the name afterward to avoid overwriting the file when you next run GeneSplitter. 

3. Reading Dot-Bracket Notation

    This section is not mandatory for using the program but is helpful for understanding the process.

    Dot-Bracket notation is used to describe the structure of RNA molecules in 2-dimensional space, which is useful since typically the 2-dimensional structure can exist independently of the 3-dimensional structure. The secondary structure of RNA is comprised of "stems" (unpaired sections) and "loops" (unpaired sections), with exceptions. 

    At a base level, the dot-bracket notation has three characters: ".", "(", and ")". This is enough to describe the majority of RNA secondary structure. Each character refers to one nucleotide, ordered the same as the sequence. The period describes unpaired nucleotides (loops), while the parentheses describe paired nucleotides (stems). The left and right parentheses describe which "side" the corresponding nucleotide is on. If the RNA structure was completely straight, the left parentheses would refer to stem nucleotides that are extending the structure away from the beginning, and the right parentheses would refer to stem nucleotides that are doubling back and 'returning' to the beginning. This means that each right parenthesis has a corresponding left parenthesis later in the dot-bracket notation.

    There are other characters to indicate more complicated structures (such as "+" for g-coupling and the use of slashes for pseudoknots), but those aren't supported in this program.

If following the instructions above, you should now be able to use the CutSequenceFinder package. 

## Reference Guide

If using this program, it is recommended to cite the following:

(This will be filled in once the paper containing this program has been published)

## License

Please read the copyright notice in the file (Insert License Here)!

## Contact

(Needs to be filled in)

[vrna_website]: https://www.tbi.univie.ac.at/RNA
[vrna_github]: https://github.com/ViennaRNA/ViennaRNA
[vrn_package]: https://www.tbi.univie.ac.at/RNA/#
[RNAfold]: http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi
[Genbank]:https://www.ncbi.nlm.nih.gov/genbank/
