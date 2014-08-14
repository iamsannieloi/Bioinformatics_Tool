import numpy as np
import matplotlib.pyplot as plt
import math
from collections import defaultdict


def findavg(pvalue):
    avg = [0,0]
    for i in range(2,len(pvalue)-2):
        total = np.mean(pvalue[i-2:i+3])
        avg.append(total)
    return avg
        
        
def plotmotif(termlist, hyper_file, overlay, cutoff, coeffvalue):
    
    #open files to be plotted
    f = open(hyper_file,'r')
    f = f.readlines()
    g = open(termlist, 'r')
    t = g.readlines()
    
    counter = 0 #used to keep track of the position
    pvalue = [ [] for i in range(len(t))]
    motif = [ [] for j in range(len(t))]
    for term in t: #loops for the terms listed that wanted to be graphed
        m = term.replace('\n','')
        
        for line in f: # loops for the motifs that are in the hypergeometric file
            a = line.replace('\n','')
            b = a.split('\t')
            c = b[1].split('|')
            
            dicterm = {}
            for var in c[0:cutoff]: # loops for the terms for that particular motif
                d = var.split(', ')
                dicterm[d[0]] = d[1]
            # I get a math domain error if the number is too low, so I set a threshold where if the number is lower than 1e-49, the pvalue defaults to 1e-49
            # If the term was not in the cutoff, then it gets set a pvalue of 1
            if m in dicterm:
                if float(dicterm[m]) >= 1e-50:
                    pvalue[counter].append(-1 *np.log10(float(dicterm[m]) * 348 * 10905 ))
                else:
                    pvalue[counter].append(-1 *np.log10(1e-49 * 348 * 10905))
                motif[counter].append(b[0])
            else:
                pvalue[counter].append(-1 * np.log10(1))
                motif[counter].append(b[0])
                
        counter += 1

    #styles of lines are created
    color = ['b','r','g','m','c','y','k'] 
    linestyles = ['-','--','-.',':']
    style = []
    for l in linestyles:
        for c in color:
            style.append(c+l)
    
    #SUBPLOT                            
    if overlay == False or overlay == 'F' or overlay == 'f' or overlay == 'false':
        
        difflst = []
        diffmotif = []
        difflst.append(pvalue[0])
        diffmotif.append(t[0])
        counter = 0
        
        new = defaultdict(list)
        new[counter].append(t[0])
        
        #clusters the results
        for j in range(1,len(pvalue)):
            close = False
            for o in range(0, len(difflst)):
                
                if np.corrcoef(pvalue[j], difflst[o])[0][1] >= coeffvalue and close == False:
        
                    new[o].append(t[j])
                    close = True              
            if close == False:
                counter += 1
                new[counter].append(t[j])
                difflst.append(pvalue[j])
                diffmotif.append(t[j])
        #creates the subplots
        for k in range(len(difflst)):
           
            plt.subplot(len(difflst),1,k)
            #c = plt.plot(difflst[k],color[k], label = new[k][0])
            c = plt.plot(difflst[k], 'k', label = new[k][0])
            plt.plot(findavg(difflst[k]),'r')
            
            for g in range(1,len(new[k])):
                c = plt.plot(difflst[k], 'w', label = new[k][g], alpha = 0)
        
            
            plt.subplots_adjust(hspace=0.5)
            plt.title(diffmotif[k].replace('\n',''), fontsize = 10)
            plt.xlabel('Motif Position')
            plt.ylabel('-log(p-value)')
            plt.ylim(ymin=0)
            plt.subplots_adjust(right = 0.6)
            plt.legend(bbox_to_anchor=(1,1), fontsize = 'x-small', ncol = 2, loc =2, borderaxespad = 0)
            
        plt.show()
    
    #OVERLAY
            
    elif overlay == True or overlay == 'T' or overlay == 't' or overlay == 'true':
        
        plt.hold(True)
        
        diff = []
        counter = 0
        plt.plot(pvalue[0],style[counter], label = t[0].replace('\n',''))
        diff.append(pvalue[0])
        for k in range(1,len(pvalue)):
            similar = False
            
            for v in diff:
            
                if np.corrcoef(pvalue[k],v)[0][1] >= coeffvalue and similar == False:
                    plt.plot(pvalue[k], color = '0.5', alpha = 0.25)
                
                    #plt.plot(pvalue[k], color[counter], label = t[k].replace('\n',''), alpha = 0.0)
                    similar = True
            if similar == False:
                counter += 1
                plt.plot(pvalue[k],style[counter], label = t[k].replace('\n',''))
                diff.append(pvalue[k])
            
        legend = plt.legend()
        legend.get_frame().set_alpha(0.5)
        plt.title('Correlated Coefficient: ' + str(coeffvalue) + ' Cutoff: ' + str(cutoff))
        plt.xlabel('Motif position')
        plt.ylabel('-log(P-value)')
        plt.ylim(ymin = 0)
        plt.show()
             
        
        
        
        
