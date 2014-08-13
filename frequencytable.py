##@author: Sannie Loi
##This program takes a text file produced from hypergeometric program and produces a text file 
## where the first column is the term, second column is frequency, and third column is all the 
## motifs associated with that term
##
##example: freqtable(TopPvalue.txt,'freqtable.txt',4)

##**If the last parameter is higher than the number of terms present, then it will default to the
##  max term 

from collections import OrderedDict
from collections import defaultdict

#this turns the list into a str that for writing to file
def turntostr(motiflist):
    motstr = motiflist[0]
    for i in motiflist[1:]:
        motstr = motstr + ', ' + i
    return motstr
        
#this counts the number of times the term occurs and has a list of the motif that is associated with
#that term
def freqtable(hypofile, outputfile, topNterm):
    
    new = open(hypofile, 'r')
    
    out = open(outputfile,'w')
    
    #Here I am creating two different dictionaries, one for frequency, one for motif
    
    count = {} #keep track of frequency(key = term, value = frequency)
    
    motif = defaultdict(list) #keep track of motif 
                              #(key = term, value = list of motif)
                              
                              
    for line in new:
        #this takes a line from the read file and splits the info
        s = line.replace('\n','')
        d = s.split('\t') #d[0] = motif, d[1] = all term for that motif
        topterm = d[1].split('|')
        
        #if topNterm > len(topterm):
        #   out.write('topNterm is too high, the highest it can go is: ' + str(len(topterm)))
            
    
        for rank in topterm[0:topNterm]:
            termpvalue = rank.split(', ')#termpvalue[0] = term, termpvalue[1] = p-value
        
            if termpvalue[0] in count: #occurs if key found in dictionary
                count[termpvalue[0]] += 1
                motif[termpvalue[0]].append(d[0])
            else: #occurs if not a key in the dictionary
                count[termpvalue[0]] = 1
                motif[termpvalue[0]].append(d[0])
    
    #sorts the dictionary by frequency (highest to lowest)
    freq = OrderedDict(sorted(count.items(), key=lambda x: x[1], reverse = True))
    
    #writes to file
    
    for i in freq:
        out.write(i + '\t' + str(freq[i]) + '\t' + turntostr(motif[i]) + '\n')
    out.close()
    new.close()