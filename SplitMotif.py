
import re

def findmotif(frag, seq, window):
    
    #this loops for as long as the sequence length          
    for i in range(0, len(seq)-window+1):
        counter = 0
        
        #this loop compares the fragment with the window from the sequence and looks for a match with an allowance of 1 mismatch        
        for j in range(0, len(frag)):
            
            if frag[j] == seq[i:i+window][j]:
                counter = counter + 1
                
        if counter >= window-1: #allows for one mismatch
            return True
            
    return False
                     
                                               
def slide(seq, window):    
    counter = 0
    mini = ""
    mainlist = []
    for i in range(0,len(seq)-window+1):
        mainlist.append(seq[i:i+window]) 
    
    return mainlist
        
        


def splitMotif(slidingwindow, userfile, outputname):
    
    file = open(userfile,'r')
    lst = file.readlines()
    term = ''
    for i in lst:
        oneline = i.replace('\n','')
        term = term + oneline
        
    
    proteome = open("TAIR10_pep_20101214.txt", 'r')
    proteome1 = proteome.readlines()
    
    
    fragmentlist = slide(term,slidingwindow)
    
    new = open(outputname,'w')
    
    #here gets the sequences from the text file 
    biglist = []
    somestr = []
    for frag in fragmentlist:
        
        motiflist = []
        
        for line in proteome1:
            if line[0] != '>':             
                somestr.append(line)
                   
            elif line[0] == '>': 
                str1 = ''.join(somestr) 
                finish = str1.replace('\n','')
             
                #if the sequence contains the motif this would occur       
                if findmotif(frag.upper(), finish, slidingwindow) == True: 
                    
                    #writes the GI accession line in the new file
                    #new.write(GI + '\n')
                    motiflist.append(GI)
                #here it refreshes somestr and assigns GI the accession num    
                somestr = [] 
                GI = line[1:10]
        if motiflist != []:
            t = ','.join(motiflist)
            new.write(frag + '\t' + t + '\n')
            
    new.close()
