
import decimal
#given a motif file with a list of arabidopsis protein 
#calculate the hypergeometric distribution using the GO_detail file
#column 1 is the terms used
#column 11 is the protein that is associated with that term

#helper function for binomial theorem
def choose(n,k):
    if 0 <= k <= n:
        nvar = 1
        kvar = 1
        for t in xrange(1,min(k,n-k)+1):
            nvar *= n
            kvar *= t
            n -= 1
        g = decimal.Decimal(nvar // kvar)
        return g
    else:
        return 0

#Population 
#Making a list of all the protein in the population and counting the total population
def hypergeometric(userfile, output):
    progenelist = []
    numpopulation = 0
    genelist = open('TAIR10_pep_20101214.txt', 'r')
    for gene in genelist:
        if gene[0] == '>':
            numpopulation += 1
            progenelist.append(gene[1:10])

    #Reads the GO file
    file = open("Ara.DetailInfo",'r')
    x = file.readlines()


    motiflist = open(userfile,'r')
    finish = open(output,'w')
    #Gets the protein list for a motif (sample) and checks number of matches in the terms
    for motif in motiflist:
        s1 = motif.replace('\n','')
        m = s1.split('\t')
        protlist = m[1].split(',')
    
        biglist = []
        #loop for EACH term
        for line in x:
       
            s = line.replace('\n','')
            something = s.split('\t')
            termprot = something[10].split(',')
    
            numsample = 0
            #loop to compare the motif with term
            #THIS IS THE NUMBER OF SUCCESSES IN THE SAMPLE
            samplematch = 0
        
            for protein in protlist:
                numsample += 1
                fix = protein.replace('\n','')
                if fix in termprot:
                    samplematch += 1
            
            #need to find number of matches in the tair file for each term
            # THIS IS NUM OF SUCCESSES IN THE POPULATION

            nummatch = 0
            for g in progenelist:
                if g in termprot:
                    nummatch += 1
    
            #CALCULATE P-VALUE
            p = (choose(nummatch, samplematch) * choose(numpopulation-nummatch, numsample - samplematch))\
                         / choose(numpopulation, numsample) 
       
            biglist.append([something[0],p])
        
        new = sorted(biglist, key = lambda pvalue: pvalue[1])
        finish.write(m[0] + '\t')
        for b in new[0:11]:
            finish.write(str(b))
        finish.write('\n')


    finish.close()

    

