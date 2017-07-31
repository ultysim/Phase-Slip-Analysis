
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as plt
import numpy as np
import scipy.optimize
import pylab
from scipy import interpolate
from scipy.fftpack import fft
import SG
import gc
from scipy.optimize import leastsq
from matplotlib.gridspec import GridSpec


#Initialzie local data storage
data = list()
prune = list()
results = list()
bestFits = list()

bval52 = [5.5845,5.5840,5.5835,5.5830,5.5825,5.5820,5.5815,5.5810,5.5805,5.5800,5.5795,5.5790,5.5785,5.5780,5.5775]
bval73 = [5.7485,5.7480,5.7475,5.7470,5.7465,5.7460,5.7455,5.7450,5.7445,5.7440,5.7435,5.7430,5.7425,5.7420,5.7415,5.7410,5.7405]
bval2011 = [5.536,5.5355,5.535,5.5345,5.534,5.5335,5.533,5.5325,5.532,5.5315,5.531,5.5305,5.53,5.5295,5.529,5.5285,5.528,5.5275,5.527,5.5265,5.526,5.5255]

cdict = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.5,0.0,0.0),
                   (1.0, 1.0, 1.0)),

         'blue':  ((0.0, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
        
yellowredblue = LinearSegmentedColormap('YellowRedBlue', cdict)

#set and get functions for data, prune, and results
def getData(i=-1):
    if i == -1:    
        return data
    else:
        return data[i]
def setData(new):
    global data
    data = list(new)
    
def getPrune(i=-1):
    if i == -1:    
        return prune
    else:
        return prune[i]  
def setPrune(new):
    global prune
    prune = list(new)
    
def getResults(i=-1):
    if i == -1:    
        return results
    else:
        return results[i]   
def setResults(new):
    global results
    results = list(new)
    
#initialzie data,prune,and results
def parseData():
    f = open('filenames')
    raw = f.read()
    parse = raw.split('\n',-1)
    out = list()
    f.close()
    for i in parse:
        if i == '':
            continue
        f = open('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\'+str(i), 'r')
        raw = f.readline()
        parse = raw.split('\r',-1)
        data =np.ndarray(shape=(2,len(parse)-2))   
        count = 0
        for a in parse: 
            if count > 0 and count < (len(parse)-1):
                hold = a.split('\t',-1)                    
                data[0,count-1]=float(hold[0])
                data[1,count-1]=float(hold[1])
            count += 1
        if data[0,0]>data[0,-1]:
            data = np.fliplr(data)
        #clean the data
        def minFunc(a):
            return max(data[1] + a*data[0]) - min(data[1] + a*data[0])
        slope = scipy.optimize.fmin(minFunc,0,disp=False)
        data[1]=data[1]+slope[0]*data[0]
        out.append(data)
    setData(out)
    
def parsePrune(filetype):
    f = open('prunenames')
    raw = f.read()
    parse = raw.split('\n',-1)
    out = list()
    f.close()
    for i in parse:
        f = open('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\'+filetype+'\\'+str(i), 'r')
        raw = f.read()
        parse = raw.split('\n',-1)
        data =np.ndarray(shape=(2,len(parse)-1))   
        count = 0
        for a in parse: 
            if count < (len(parse)-1):                               
                data[0,count]=float(a.split(' ',-1)[0])
                data[1,count]=float(a.split(' ',-1)[1])
            count += 1    
        if data[0,0]>data[0,-1]:
            data = np.fliplr(data) 
        out.append(data)
    setPrune(out)
#Creates array of all data sets with given results
#Indexed as [file][fit array][step in fit array]
def parseResults():
    f = open('resultsnames')
    raw = f.read()
    parse = raw.split('\n',-1)
    out = list()
    f.close()
    for i in parse:
        x = np.load(str(i))        
        out.append(x)
    finalhold = list()
    for i in range(len(out[0])):
        datahold = list()        
        for j in range(len(out)):
            datahold.append(out[j][i])
        finalhold.append(datahold)
    setResults(finalhold)        

#Creates a npy array file where the indeces are always scaled to original data        
def makeResults(name,unpruned=True,pruneType='soft'):
    mat = scipy.io.loadmat('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\output.mat')
    results = mat['result']    
    if unpruned == True:
        for i in results:
            for j in i:
                j[7] -= 1
                j[8] -= 1
                j[10] -= 1
        np.save(name,results)
    else:
        index = 0
        parseData()
        parsePrune(pruneType)        
        for i in results:
            data = getData(index)
            pruned = getPrune(index)             
            for j in i:
                j[7] = list(data[0]).index(pruned[0][j[7]-1])
                j[8] = list(data[0]).index(pruned[0][j[8]-1])
                j[10] = list(data[0]).index(pruned[0][j[10]-1])
            index += 1
        np.save(name,results) 
#Generates slip widths for all results files 
#Returns array as large as results, where each set of results has an array for each file 
def slipWidth(n=3):
    results = getResults()
    out = list()    
    for i in results:       
        filehold = list()
        for j in i:                   
            slips = j[:,9]
            resulthold = list()
            start = 0
            stop = 0
            sweeping = False
            i= 0
            while(i<119):
                if (not sweeping and (slips[i]>0.20 and slips[i]<1.80)):          
                    start = i
                    sweeping = True
                if sweeping and ((np.abs(slips[i]-slips[i+1])>0.02) or i==118 or not (slips[i+1]>0.20 and slips[i+1]<1.80)):
                    stop = i
                    if stop - start > n:
                        resulthold.append([start,stop+1]) #store results for the individual files
                    sweeping = False
                    i = i + 1
                else:
                    i = i + 1
            filehold.append(resulthold)#combine fileholds for the whole result
        out.append(filehold)#combine all the resultholds            
    return out
    
#finds slip with greatest fit and physical params, and returns array per file with all applicable fits          
def bestFit(n=3):  
    widths = slipWidth(n)    
    out = list()    
    indexi = 0    
    for i in widths: #Goes through the data sets
        filehold = list()
        indexj = 0        
        for j in i: #Goes through the fit results            
            for k in j: #Goes through the individual widths
                singlewidth = getResults()[indexi][indexj][k[0]:k[1]] #selects files between [start,stop+1)
                singlewidth = sorted(singlewidth, key = lambda f:f[6])
                amp = 1.5*np.std(prune[indexi][1])
                #ampstd = np.std(clip[:,1])
                ampstd = amp*.20             
                height = np.mean(prune[indexi][1])
                #heightstd = np.std(clip[:,4])
                heightstd = amp*.10             
                freq = 0.6826
                freqstd = 0.13
                while len(singlewidth) > 0:                                     
                    temp = singlewidth.pop()
                    noiseclip = (temp[9]>0.1 and temp[9]<1.9)
                    heightbool = np.abs(temp[4]-height)<heightstd 
                    ampbool = np.abs(temp[1]-amp)<ampstd
                    freqbool = np.abs(temp[0]-freq)<freqstd
                    if heightbool and ampbool and freqbool and noiseclip:
                        filehold.append(temp)
                        break
            indexj += 1
        out.append(np.array(filehold))
        indexi += 1    
    return out

#Takes index values from fits and computes original voltage values                    
def correctedFit(n=3):
    global bestFits
    data = getData()
    out = bestFit(n)
    fileindex = 0
    for i in out: #iterate through data files                
        for j in i:
            try:
                j[8] = data[fileindex][0][j[8]]
            except:
                j[8] = data[fileindex][0][j[8]-1]
            j[7] = data[fileindex][0][j[7]]
        fileindex += 1
    bestFits = out

#Returns left and right fit values for a voltage x
def yValue(fit,x):
    xa = fit    
    freq = xa[0]
    amp = xa[1] 
    phasel = xa[2]
    phaser = xa[3]
    height = xa[4]
    slope = xa[5]
    v0 =  xa[11]
    yl = height + amp*np.sin(freq*x + phasel) + slope*(x - v0)
    yr = height + amp*np.sin(freq*x + phaser)  + slope*(x - v0)        
    return yl,yr

#Takes an x axis index value and computes applicaple fit values with associted slip
#Returns an array for individual files where [ylefti,yrighti,ylefti+1,yrighti+1,dphi]
def fitValue(filenum,index,pruned=False):       
    fits = bestFits[filenum] #Want voltage values so interchangeable with data and prune
    datahold = list()
    if pruned:
        datahold = getPrune(filenum)
    else:
        datahold = getData(filenum)    
    out = list()
    x0 = datahold[0][index]
    x1 = datahold[0][index+1]
    for i in fits:
        if x0 >= i[7] and x0 <= i[8]:
            yl0,yr0 = yValue(i,x0)
            yl1,yr1 = yValue(i,x1)            
            out.append([yl0,yr0,yl1,yr1,i[9]])    
    return out 

def closestFit(filenum,index,pruned=False,mindiff=0.02):
    fits = fitValue(filenum,index,pruned) #get y values for fits that are applicable for x value
    datahold = list()
    if pruned:
        datahold = getPrune(filenum)
    else:
        datahold = getData(filenum)
    y0 = datahold[1][index]
    y1 = datahold[1][index+1]
    dy = mindiff
    hold = list()
    for i in fits:
        if np.abs(y0-i[0]) <= dy and np.abs(y1-i[3]) <= dy:
            hold.append([np.abs(y0-i[0])+np.abs(y1-i[3]),i[4]])
        elif np.abs(y0-i[1]) <= dy and np.abs(y1-i[2]) <= dy:
            hold.append([np.abs(y0-i[1])+np.abs(y1-i[2]),i[4]])
    sorthold = sorted(hold, key = lambda f:f[0])
    if len(sorthold)>0:
        return sorthold[0][1]

#Creates histogram from given data file number, allows prune to be on and change of fit thickness    
def makeHistogram(filenum,pruned=False,minv=1,maxv=1,mindiff=0.005):
    histogram = list()
    datahold = list()
    if pruned:
        datahold = getPrune(filenum)
    else:
        datahold = getData(filenum)
    if minv > 0:
        minv = datahold[0,0] + 0.05
    if maxv > 0:
        maxv = datahold[0,-1] - 0.05
    for i in range(len(datahold[0])-1):
        if datahold[0][i] > minv and datahold[0][i] <maxv: #Remove the edges or open desired window
            pass
        else:
            continue
        histogram.append(closestFit(filenum,i,pruned,mindiff))
    return histogram

#Plots histogram for given file number and histogram output 
def plotHistogram(filenum,histogram,bins=.03,clr='blue'):        
    if len(histogram)>0:    
        bins = np.arange(0,2.0,bins)
        plt.pyplot.figure(filenum)
        plt.pyplot.hist(histogram,bins,color=clr)


#Initialize all data sources        
def initialize(n=3,pruneType='soft'):
    parseData()
    parseResults()
    parsePrune(pruneType)    
    correctedFit(n)
    
#Plot a fit for a given result number:
def plotfit(filenum,resultnum,step,back=True,max=False):
    xa = getResults(filenum)[resultnum][step]
    data = getData(filenum)
    freq = xa[0]
    amp = xa[1] 
    phasel = xa[2]
    phaser = xa[3]
    height = xa[4]
    slope = xa[5]
    posleft = xa[7] 
    posright = xa[8] 
    pos =  xa[10] 
    v0 = xa[11]
    xposl = list()
    xposr = list()       
    xposl = data[0][posleft:pos+1]
    xposr = data[0][pos:posright+1]
    if max == True:
        xposl = data[0][0:pos+1]
        xposr = data[0][pos:]
    if back ==True:
        pylab.plot(data[0],data[1])
    ytargll = height + amp*np.sin(freq*xposl + phasel) + slope*(xposl - v0)
    ytarglr = height + amp*np.sin(freq*xposr + phasel) + slope*(xposr - v0)
    ytargrl = height + amp*np.sin(freq*xposl + phaser) + slope*(xposl - v0)
    ytargrr = height + amp*np.sin(freq*xposr + phaser) + slope*(xposr - v0)   
    pylab.plot(xposl,ytargll,'ko',markersize=4,mew=0)
    pylab.plot(xposr,ytarglr,'ko',markersize=1,mew=0)
    pylab.plot(xposr,ytargrr,'ro',markersize=4,mew=0)
    pylab.plot(xposl,ytargrl,'ro',markersize=1,mew=0)
#takes a result array as fitbit
def plotbit(filenum,fitbit,back=True,max=False):
    xa = fitbit
    data = getData(filenum)
    freq = xa[0]
    amp = xa[1] 
    phasel = xa[2]
    phaser = xa[3]
    height = xa[4]
    slope = xa[5]
    posleft = xa[7] 
    posright = xa[8] 
    pos =  xa[10] 
    v0 = xa[11]
    xposl = list()
    xposr = list()       
    xposl = data[0][posleft:pos+1]
    xposr = data[0][pos:posright+1]
    if max == True:
        xposl = data[0][0:pos+1]
        xposr = data[0][pos:]
    if back ==True:
        pylab.plot(data[0],data[1])
    ytargll = height + amp*np.sin(freq*xposl + phasel) + slope*(xposl - v0)
    ytarglr = height + amp*np.sin(freq*xposr + phasel) + slope*(xposr - v0)
    ytargrl = height + amp*np.sin(freq*xposl + phaser) + slope*(xposl - v0)
    ytargrr = height + amp*np.sin(freq*xposr + phaser) + slope*(xposr - v0)   
    pylab.plot(xposl,ytargll,'ko',markersize=4,mew=0)
    pylab.plot(xposr,ytarglr,'ko',markersize=1,mew=0)
    pylab.plot(xposr,ytargrr,'ro',markersize=4,mew=0)
    pylab.plot(xposl,ytargrl,'ro',markersize=1,mew=0)
    
#If there are errors with index out of bound with results, removes the 1 offest from matlab    
def fixResults():
    f = open('resultsnames')
    raw = f.read()
    parse = raw.split('\n',-1)
    f.close()
    for i in parse:       
        x = np.load(str(i))
        if x[0,0,7]>0:
            for j in x:
                j[7] = j[7] - 1
                j[8] = j[8] - 1
                j[11] = j[11] - 1
            name = str(i)+'new'
            np.save(name,x)
            
def slipAverage():
    hold = []
    out = []
    res = getResults()
    for i in res:
        for j in i:
            hold.append(j[:,9])            
        filedatalist = []
        for l in range(120):
            temp = []
            for k in range(len(hold)):
                if hold[k][l] > .1 and hold[k][l] < 1.9:
                    temp.append(hold[k][l])
                else:
                    temp.append(0)
            filedatalist.append(temp)
        means = []        
        for l in filedatalist:
            mean = np.mean(l)
            std = np.std(l)            
            sums = 0
            count = 0            
            for k in l:
                meanbool = np.abs(k-mean)<std
                if meanbool:
                    sums += k
                    count +=1            
            if count > 0:
                means.append(sums/count)
            else:
                means.append(0)
        out.append(means)
    return out
            
def smoothSlip():
    results = list(getResults())
    out = list()
    n = .02    
    for i in results:       
        filehold = list()
        for j in i:                   
            slips = j[:,9]
            avg = list()
            start = 0
            stop = 0
            sweeping = False
            i= 0
            while(i<119):
                if not sweeping and (slips[i]>0.20 and slips[i]<1.80) and np.abs(slips[i]-slips[i+1])<n:          
                    start = i
                    sweeping = True
                    avg.append(slips[i])
                    continue
                if sweeping and np.abs(slips[i]-np.mean(avg))<n and (i != 118 or (slips[i]>0.20 and slips[i]<1.80)):
                    avg.append(slips[i])
                    i += 1
                elif np.abs(slips[i-1]-np.mean(avg))<n and np.abs(slips[i+1]-np.mean(avg))<n and (i != 118 or (slips[i]>0.20 and slips[i]<1.80)) and np.abs(slips[i]-np.mean(avg))<.1:
                    i += 1
                else:
                    stop = i
                    sweeping = False
                    if len(avg)>0:
                        slips[start:stop] = np.mean(avg)
                    i += 1
                    avg = list()
            filehold.append(slips)#combine fileholds for the whole result
        out.append(filehold)#combine all the resultholds            
    return out

def modifySlips():
    slips = smoothSlip()
    results = list(getResults())
    indexi = 0
    for i in results:
        indexj = 0
        for j in i:                   
            start = 0
            stop = 0
            sweeping = False
            i= 0
            while(i<119):
                if (not sweeping and (slips[indexi][indexj][i]>0.20 and slips[indexi][indexj][i]<1.80)):          
                    start = i
                    sweeping = True
                if sweeping and ((np.abs(slips[indexi][indexj][i]-slips[indexi][indexj][i+1])>0.02) or i==118 or not (slips[indexi][indexj][i+1]>0.20 and slips[indexi][indexj][i+1]<1.80)):
                    stop = i
                    if stop - start > 2:
                        j[start:stop,7] = j[start,7]
                        j[start:stop,8] = j[stop-1,8]
                        j[start:stop,9] = j[start,9]
                    sweeping = False
                    i = i + 1
                else:
                    i = i + 1
            indexj += 1
        indexi += 1
    return results
    

#Uses average chi square for left and right fits and returns best fit for small range    
def qualFit(filenum,index,width=0.5,pruned=True):
    fit = bestFits[filenum] #Want voltage values so interchangeable with data and prune
    data = list()
    if pruned:
        data = getPrune(filenum)
    else:
        data = getData(filenum)     
    xlow = search(data,data[0][index]-width)
    xhigh = search(data,data[0][index]+width)
    bfhold = list()
    if len(fit)==0:
        return None
    for a in fit:
        if data[0][index] >= a[7] and data[0][index] <= a[8]:
            bfhold.append(a)
    if len(bfhold)==0:
        return None
    qual = list()
    for i in fit:
        out = list()
        for j in range(xlow,xhigh+1):
            out.append((data[1][j]-yValue(i,data[0][j]))**2/yValue(i,data[0][j]))
        final = sum(out)
        qual.append(final[0]+final[1])
    return fit[qual.index(min(qual))]
#Checks if point is on a fit and returns the delta phi and left, right, or mid
#returns [either delta phi or -1,0 for nothing 1 left 2 right,x value, and fit quality] 
def onFit(filenum,fitnum,index,pruned=False,mindiff=0.030):
    i = index
    datahold = list()
    if pruned:
        datahold = getPrune(filenum)
    else:
        datahold = getData(filenum)    
    x = datahold[0][i]
    y = datahold[1][i]
    if len(bestFits[filenum])==0:
        return [-1,0,x,0]    
    fit = bestFits[filenum][fitnum] #Want voltage values so interchangeable with data and prune
    dy = mindiff
    if x >= fit[7] and x <= fit[8]:    
        left = np.abs(y-yValue(fit,x)[0])<dy
        right = np.abs(y-yValue(fit,x)[1])<dy
        if np.abs(yValue(fit,x)[0]-yValue(fit,x)[1])<2*dy:
            return [-1,0,x,fit[6]]
        if left == False and right == False:
            return [fit[9],0,x,fit[6]]
        if left:
            return [fit[9],1,x,fit[6]]
        if right:
            return [fit[9],2,x,fit[6]]
    return [-1,0,x,0]

def onFitSweep(filenum,pruned=True,minv=1,maxv=1,mindiff=0.030):
    datahold = list()
    out = list()
    if pruned:
        datahold = getPrune(filenum)
    else:
        datahold = getData(filenum)
    if minv > 0:
        minv = datahold[0,0] + 0.05
    if maxv > 0:
        maxv = datahold[0,-1] - 0.05
    for j in range(len(bestFits[filenum])):
        hold = list()
        i = 0
        while i < (len(datahold[0])-1):
            if datahold[0][i] > minv and datahold[0][i] <maxv: #Remove the edges or open desired window
                pass
            else:
                i += 1 
                continue
            hold.append(onFit(filenum,j,i,pruned,mindiff))
            i += 1
        out.append(np.array(hold))
    return np.array(out)
            

def smoothOnFit(filenum,pruned=True,minv=1,maxv=1,mindiff=0.030):
    ofs = np.array(onFitSweep(filenum,pruned,minv,maxv,mindiff))
    for j in ofs:
        i = 0
        while i < len(j):
            if j[i][1] == 1 or j[i][1] == 2:
                start = i
                state = j[i][1]
                i += 1
                while True and i < len(j):
                    if j[i][0] < 0:
                        i += 1                    
                        break
                    if j[i][1] == state:
                        if np.abs(j[start][2]-j[i][2])<0.04:
                            j[start:i,1] = state
                        break
                    elif j[i][1] == 0:
                        i += 1
                    else:
                        i += 1
                        break
            else:
                i += 1
    return ofs
    

#Finds edges of data for fit, returns dphi,left edge, right edge, dV of edges
def findEdge(filenum,dV=0.025,pruned=True,minv=1,maxv=1,mindiff=0.030):                
    sof = smoothOnFit(filenum,pruned,minv,maxv,mindiff)    
    out = list()
    for j in sof:
        i = 0
        hold = list()
        while i < len(j):
            state = j[i,1]
            if state == 1 or state == 2:
                start = i
                while True:
                    i += 1
                    if i >= len(j):
                        break
                    if j[i,1] == state:
                        continue
                    elif j[i,0] < 0:
                        break                        
                    elif np.abs(j[start,2]-j[i-1,2]) > dV:
                        hold.append([j[start,0],j[start,1],j[start,2],j[i-1,2],np.abs(j[start,2]-j[i-1,2]),j[start,3]])
                        break
                    else:
                        break
            else:
                i += 1
        out.append(hold)
    return out
#If a slip exists, returns the right point of edge and left point of next edge
def findSlips(filenum,dV=0.025,pruned=True,minv=1,maxv=1,mindiff=0.030,minsep=0.2):
    fe = findEdge(filenum,dV,pruned,minv,maxv,mindiff)
    out = list()
    for k in fe:
        hold = list()
        for i in range(len(k)-1):
            j = i + 1
            if k[i][1] != k[j][1]:
                if np.abs(k[i][3]-k[j][2])<minsep:
                    hold.append([k[i][0],k[i][1],k[i][3],k[j][2],k[i][5]])
        out.append(np.array(hold))
    return np.array(out)

def finalSlips(filenum,dV=0.025,pruned=True,minv=1,maxv=1,mindiff=0.030,minsep=0.2):
    fs = findSlips(filenum,dV,pruned,minv,maxv,mindiff,minsep)
    alll = list()
    for i in fs:
        for j in i:
            alll.append(j)
    i = 0
    while i < len(alll):
        delete = False    
        j = i+1
        edge = alll[i][2]
        while j < len(alll):
            if np.abs(edge - alll[j][2]) < dV:
                if alll[j][4] <= alll[i][4]:
                    del alll[j]
                else:
                    del alll[i]
                    delete = True
                    j = len(alll)
            else:
                j += 1
        if not delete:
            i += 1
    return sortFinal(filenum,alll)

def sortFinal(filenum,fSOutput):
    alll = list(fSOutput)    
    out = list()
    bf = bestFits[filenum]
    for j in bf:
        hold = list()
        fit = j[6]
        i = 0
        while i <len(alll):
            if alll[i][4] == fit:
                hold.append(alll[i])
                del alll[i]
            else:
                i +=1
        out.append(np.array(hold))
    return np.array(out)
    
def countSlips(filenum,dV=0.025,pruned=True,minv=1,maxv=1,mindiff=0.030,minsep=0.2):
    fs = finalSlips(filenum,dV,pruned,minv,maxv,mindiff,minsep)
    out = list()
    for i in fs:
        for j in i:
            out.append(j[0])
    return out

#Stores the location of phase slip events, takes finalSlip and returns a list of voltage locations
#and phase values.
def countEdge(filenum,dtype='52',dV=0.025,pruned=True,minv=1,maxv=1,mindiff=0.030,minsep=0.2):
    fs = finalSlips(filenum,dV,pruned,minv,maxv,mindiff,minsep)
    hold = list()
    if dtype=='52':
        b = bval52
    elif dtype=='73':
        b = bval73
    else:
        b = bval2011
    for i in fs:
        for j in i:
            if np.abs(j[2]-j[3])<=2:
                hold.append([(j[2]+j[3])/2,b[filenum],j[0]])
            else:
                hold.append([j[2]+0.05,b[filenum],j[0]])
    return np.array(hold)

#Sorts the edges into boxes defined by the user
def sortEdge(boxes,dtype='52',dV=0.025,pruned=True,minv=1,maxv=1,mindiff=0.030,minsep=0.2):
    hold = list()
    final = list() 
    for i in range(len(prune)):
        hold.extend(countEdge(i,dtype,dV,pruned,minv,maxv,mindiff,minsep))
    for i in boxes:
        store = list()
        j = 0
        while j < len(hold):
            if hold[j][2] >= i[0] and hold[j][2] <= i[1]:
                store.append(np.array(hold.pop(j)))
            else:
                j+=1        
        final.append(np.array(store))
    final.append(np.array(hold))
    return final        

    
    
#Runs entire fit algorithm, computes histograms for each data file and a full histogram    
def mainRun(n=3,bins=0.03,dV=0.025,pruned=True,pruneType='soft',minv=1,maxv=1,mindiff=0.030,minsep=0.2):
    datahold = list()
    full = list()
    if pruned:
        datahold = getPrune()
    else:
        datahold = getData()    
    for i in range(len(datahold)):
        hist = countSlips(i,dV,pruned,minv,maxv,mindiff,minsep)
        full.append(hist)
        plotHistogram(i,hist,bins)
    out = [item for sublist in full for item in sublist]
    plotHistogram('full',out,bins)
    return out

def selectRun(files,n=3,bins=0.03,dV=0.025,pruned=True,pruneType='soft',minv=1,maxv=1,mindiff=0.030,minsep=0.2):
    datahold = list()
    full = list()
    if pruned:
        datahold = getPrune()
    else:
        datahold = getData()    
    for i in range(len(datahold)):
        if i in files:
            hist = countSlips(i,dV,pruned,minv,maxv,mindiff,minsep)
            full.append(hist)
    out = [item for sublist in full for item in sublist]
    plotHistogram('full',out,bins)
    return out

    
def plotEdges(filenum,edges):
    fs = edges
    data = getData(filenum)
    for i in range(len(fs)):
        if len(fs[i]) == 0:
            continue
        pylab.figure(i)
        plotPrune(filenum,bestFits[filenum][i])
        xpoints = list()
        for j in fs[i]:
            xpoints.append(j[2])
            xpoints.append(j[3])
        ypoints = list()
        for j in xpoints:
            ypoints.append(data[1][list(data[0]).index(j)])
        pylab.scatter(xpoints,ypoints,color='black',s=300)
        pylab.show()
        
def teleLength(filenum,dV=0.025,pruned=True,minv=1,maxv=1,mindiff=0.030,minsep=0.2):
    fs = finalSlips(filenum,dV,pruned,minv,maxv,mindiff,minsep)
    out = list()
    for i in fs:
        for j in range(len(i)-1):
            if np.abs(i[j][3]-i[j+1][2]) < 1:
                out.append(np.abs(i[j][3]-i[j+1][2]))
    return out
 
def teleCount(filenum,index,pruned=False,mindiff=0.025):
    i = index    
    fits = bestFits[filenum] #Want voltage values so interchangeable with data and prune
    datahold = list()
    if pruned:
        datahold = getPrune(filenum)
    else:
        datahold = getData(filenum)    
    x = datahold[0]
    y = datahold[1]
    dy = mindiff
    bfhold = list()
    if len(fits)==0:
        return i+1,None
    for a in fits:
        if x[i] >= a[7] and x[i] <= a[8]:
            bfhold.append(a)
    if len(bfhold)==0:
        return i+1,None
    bfhold = sorted(bfhold, key = lambda f:f[6])
    best = bfhold.pop()    
    left = np.abs(y[i]-yValue(best,x[i])[0])<dy
    right = np.abs(y[i]-yValue(best,x[i])[1])<dy
    if np.abs(yValue(best,x[i])[0]-yValue(best,x[i])[1])<2*dy:
        return i+1,None
    if left == right:
        return i+1,None
    elif left:
        if np.abs(y[i]-yValue(best,x[i])[0])<dy and np.abs(y[i-1]-yValue(best,x[i-1])[0])<dy:           
            if np.abs(y[i+1]-yValue(best,x[i+1])[0])<dy:
                return i+1,None
            else:
                while True:
                    i += 1
                    if i >= len(x)-2:
                        return i, None
                    if np.abs(y[i]-yValue(best,x[i])[0])<dy:
                        return i,None
                    if np.abs(y[i]-yValue(best,x[i])[1])<dy:
                        if np.abs(np.mean([y[i],y[i+1],y[i+2]])-np.mean([yValue(best,x[i])[1],yValue(best,x[i+1])[1],yValue(best,x[i+2])[1]]))<dy:                         
                            return i,best[9]
                        else:
                            return i,None
    elif right:
        if np.abs(y[i]-yValue(best,x[i])[1])<dy and np.abs(y[i-1]-yValue(best,x[i-1])[1])<dy:           
            if np.abs(y[i+1]-yValue(best,x[i+1])[1])<dy:
                return i+1,None
            else:               
                while True:
                    i += 1
                    if i >= len(x)-2:
                        return i, None
                    if np.abs(y[i]-yValue(best,x[i])[1])<dy:
                        return i,None
                    if np.abs(y[i]-yValue(best,x[i])[0])<dy:
                        if np.abs(np.mean([y[i],y[i+1],y[i+2]])-np.mean([yValue(best,x[i])[0],yValue(best,x[i+1])[0],yValue(best,x[i+2])[0]]))<dy:                       
                            return i,best[9]
                        else:
                            return i,None                    
    return i+1,None


def teleSweep(filenum,pruned=False,minv=1,maxv=1,mindiff=0.025):
    histogram = list()
    datahold = list()
    out = list()
    if pruned:
        datahold = getPrune(filenum)
    else:
        datahold = getData(filenum)
    if minv > 0:
        minv = datahold[0,0] + 0.05
    if maxv > 0:
        maxv = datahold[0,-1] - 0.05
    i = 0
    while i < (len(datahold[0])-1):
        if datahold[0][i] > minv and datahold[0][i] <maxv: #Remove the edges or open desired window
            pass
        else:
            i += 1 
            continue
        i,out = teleCount(filenum,i,pruned,mindiff)
        if out is not None:
            histogram.append(out)        
    return histogram   
    
   

def search(data, t):
    seq = data[0]    
    min = 0
    max = len(seq) - 1
    while True:        
        m = (min + max) // 2
        if max < min:
            if np.abs(t-seq[m]) < np.abs(t-seq[m-1]):          
                return m
            else:
                return m-1        
        if seq[m] < t:
            min = m + 1
        elif seq[m] > t:
            max = m - 1
        else:
            return m
            

def clip(data,min,max):
    li = search(data,min)
    ri = search(data,max)
    return data[:,li:ri+1]
    
def plotSingle(filenum,fitbit,left,right,color,phi,style='-',size=2,back=True,max=False,prune=False):
    xa = fitbit
    data = getData(filenum)
    if prune == True:
        pruned = getPrune(filenum)
    freq = xa[0]
    amp = xa[1]
    if phi == 'l':
        phase = xa[2]
    else:
        phase = xa[3]
    height = xa[4]
    slope = xa[5]  
    v0 = xa[11]     
    xpos = clip(data,left,right)[0]
    if max == True:
        xpos = data[0]
    if back ==True:
        if prune == False:
            pylab.plot(data[0],data[1])
        else:
            pylab.plot(pruned[0],pruned[1])
    ytarg = height + amp*np.sin(freq*xpos + phase) + slope*(xpos - v0)
    pylab.plot(xpos,ytarg,color,linewidth=size,ls=style,mew=0)
    
#Takes 1d array of values and returns normalized 1d array
def normalize(data):
    data = smooth(data)
    minv = min(data)
    maxv = max(data)
    hold = np.array(data)
    for i in range(len(hold)):
        hold[i] = (hold[i]-minv)/(maxv-minv)
    return hold

def pruneRoutine(data,points=30,count=20,tol=0.072):
    new = np.ndarray(shape=(data.shape))
    new[1] = 1
    new[0] = data[0]
    j = points   
    while j < (len(new[1])-points):
        counts = 0        
        for k in range(1,points+1):
            if np.abs(data[1,j]-data[1,j+k]) > tol:
                counts += 1
            if np.abs(data[1,j]-data[1,j-k]) > tol:
                counts += 1
        if counts > count:
            new[1,j] = 0
        j += 1
    return new

#Create two column matrix for exporting
def twoCol(data):
    out = [[0] * 2 for i in range(len(data[0]))]
    for i in range(len(data[0])):
        out[i][0] = data[0][i]
        out[i][1] = data[1][i]
    return out
    
def prunedData(data,points=30,count=20,tol=0.072):
    pruned = pruneRoutine(data,points,count,tol)
    out = list(data)
    out[0] = list(out[0])
    out[1] = list(out[1])
    j = 0
    for i in pruned[1]:
        if i==0:
            del out[0][j]
            del out[1][j]
        else:
            j += 1
    return out

def createPruned(points=30,count=20,tol=0.072,name='x'):
    data = getData()
    f = open('filenames')
    raw = f.read()
    parse = raw.split('\n',-1)
    f.close()
    index = 0
    for i in parse:
        if i == '':
            continue
        pruned = prunedData(data[index],points,count,tol)
        out = twoCol(pruned)
        np.savetxt(str(i)+name,out)
        index += 1        

def interp2D(bfieldval,smoothval=501,colortype=yellowredblue,n=4096):
    parseData()   
    data = getData()
    cmlabel = 'Normalized $\Delta$R$_D$ (k$\Omega$)'
    data = [x for (y,x) in sorted(zip(bfieldval,data), key=lambda pair: pair[0])]
    bval = sorted(bfieldval)    
    xmin = None
    xmax = 0
    for i in data:
        if i[0][0] > xmin:
            xmin = i[0][0]
        if i[0][-1] < xmax:
            xmax = i[0][-1]
    xmin += 0.5
    xmax -= 0.5
    xv = np.linspace(xmin,xmax,n)
    intv = list()
    if smoothval != 0:
        for i in data:
            dataclip = clip(i,xmin,xmax)
            #dataclip[1] = SG.savitzky_golay(dataclip[1],1751,3)
            tck = interpolate.splrep(dataclip[0], dataclip[1])
            y = interpolate.splev(xv,tck)
            y = SG.savitzky_golay(y,smoothval,3)
            intv.append(normalize(y))
    else:
        for i in data:
            dataclip = clip(i,xmin,xmax)
            tck = interpolate.splrep(dataclip[0], dataclip[1])
            y = interpolate.splev(xv,tck)
            intv.append(normalize(y))
    xb = np.linspace(bval[0],bval[-1],n)
    b = list()
    for i in range(len(intv[0])):
        hold = list()
        for j in range(len(intv)):
            hold.append(intv[j][i])
        b.append(hold)
    intb = list()
    for i in b:
        tck = interpolate.splrep(bval,i)
        y = interpolate.splev(xb,tck)
        intb.append(y)
    Z = np.array(intb)
    fig = pylab.figure('ColorMap')
    X,Y = pylab.meshgrid(xb,xv)
    pylab.pcolormesh(X,Y,Z,vmin=0,vmax=1,cmap=colortype)
    pylab.ylim(xmin,xmax)
    pylab.xlim(bval[0],bval[-1])
    pylab.ylabel('Gate Voltage (mV)')
    pylab.xlabel('B Field (T)')
    cbar = pylab.colorbar()
    cbar.set_label(cmlabel)
    plt.ticker.ScalarFormatter(useOffset=None)
    pylab.show()
    gc.collect()    
    return X,Y,Z

#Smooths yvalues    
def smooth(data,count=10,tol=0.072):
    yval = list(data)
    for i in range(len(yval)):
        left = True
        right = True
        meanl = 0
        meanr = 0
        if i - count >= 0:
            meanl = np.mean(yval[i-count:i+1])
            left = np.abs(meanl-yval[i]) < tol
        if not left:
            if i + count < len(yval):
                meanr = np.mean(yval[i:i+count+1])
                right = np.abs(meanr-yval[i]) < tol
        if not left and not right:
            yval[i] = np.mean([meanl,meanr])
        if not left:
            yval[i] = meanl
        if not right:
            yval[i] = meanr
    return np.array(yval)




def noiseProfile(data,count=10,tol=0.1):
    yval = list(data)
    out = list()
    for i in range(len(yval)):
        if i - count >= 0:
            maxv = max(yval[i-count:i+1])
            minv = min(yval[i-count:i+1])
            if maxv-minv < tol:
                out.append(maxv-minv)
        elif i + count < len(yval):
            maxv = max(yval[i:i+count+1])
            minv = min(yval[i:i+count+1])
            if maxv-minv < tol:
                out.append(maxv-minv)
    return np.mean(out)

def stretch(x):
    d = x - 0.5
    f = 0.20*np.sign(x)*np.exp(-x**2/0.5)
    return d + f + 0.5
    
def plotPrune(filenum,fitbit,back=True,max=False):
    xa = fitbit
    prune = getPrune(filenum)
    freq = xa[0]
    amp = xa[1] 
    phasel = xa[2]
    phaser = xa[3]
    height = xa[4]
    slope = xa[5]
    posleft = xa[7] 
    posright = xa[8] 
    v0 = xa[11]
    xposl = list()
    xposr = list()       
    xposl = np.linspace(posleft,v0,10000)
    xposr = np.linspace(v0,posright,10000)
    if max == True:
        xposl = np.linspace(prune[0][0],v0,10000)
        xposr = np.linspace(v0,prune[0][-1],10000)
    if back ==True:
        pylab.plot(prune[0],prune[1])
    ytargll = height + amp*np.sin(freq*xposl + phasel) + slope*(xposl - v0)
    ytarglr = height + amp*np.sin(freq*xposr + phasel) + slope*(xposr - v0)
    ytargrl = height + amp*np.sin(freq*xposl + phaser) + slope*(xposl - v0)
    ytargrr = height + amp*np.sin(freq*xposr + phaser) + slope*(xposr - v0)   
    pylab.plot(xposl,ytargll,'ko',markersize=4,mew=0)
    pylab.plot(xposr,ytarglr,'ko',markersize=1,mew=0)
    pylab.plot(xposr,ytargrr,'ro',markersize=4,mew=0)
    pylab.plot(xposl,ytargrl,'ro',markersize=1,mew=0) 
    
def plotFig(data,fitbit):
    xa = fitbit
    freq = xa[0]
    amp = xa[1] 
    phasel = xa[2]
    phaser = xa[3]
    height = xa[4]
    slope = xa[5]
    posleft = xa[7] 
    posright = xa[8] 
    v0 = xa[11]     
    xpos = np.linspace(posleft,posright,10000)
    ytargl = height + amp*np.sin(freq*xpos + phasel) + slope*(xpos - v0)
    ytargr = height + amp*np.sin(freq*xpos + phaser) + slope*(xpos - v0)
    pylab.plot(data[0],data[1])
    pylab.plot(xpos,ytargl,'k',linewidth=2)
    pylab.plot(xpos,ytargr,'r',linewidth=2) 
        
#Creates subset of histogram between high and low
def histWindow(hist,low,high):
    out = list()    
    for i in hist:
        if i > low and i < high:
            out.append(i)
    return out

#Creates x,y array using midpoints of bins    
def histPoints(n,bins):
    x = list()
    y = list()
    for i in range(len(n)):
        if n[i] > 0:
            y.append(n[i])
            x.append((bins[i]+bins[i+1])/2)
    return [x,y]

    
def histFit(hist,low,high,start,stop,bins=0.03):
    hw = histWindow(hist,low,high)
    n,bins,patches = pylab.hist(hw, np.arange(low,high,bins))
    hp = histPoints(n,bins)
    xstart = float(start)
    xstop = float(stop)
    xaxis = np.linspace(xstart,xstop,1000)
    def f(x, a, b, c):
        return a * np.exp(-(x - b)**2.0 / (2 * c**2))
    popt, pcov = scipy.optimize.curve_fit(f, hp[0], hp[1])
    return xaxis,f(xaxis,*popt),popt


def primeQuart(data,xL,xH):
    xmin = search(data,xL)
    xmax = search(data,xH)    
    new = np.array([data[0][xmin:xmax],data[1][xmin:xmax]])
   
    h = sum(new[1])/len(new[1])
    a = np.std(new[1])
    f = 2*np.pi/10
    val = [h,a,f,0]


    return [new,val]
    
def firstFit(data,val,n=150):
    def f(vars,data,weight):
        x = data[0]
        yaxis = data[1]       
        h = vars[0]        
        a = vars[1]
        f = vars[2]
        p = vars[3]        
        model = (h+a*np.sin(f*x+p*np.pi))
        return (yaxis-model)*weight
    
    w = 1-np.zeros(len(data[0]))
    i = n
    for i in range(len(data[0])-n):
        s = sum(data[1][i-n:i+n])/(2*n)
        if np.abs(data[1][i]-s)>0.05:
            w[i] = 0
    
    
    result = leastsq(f, val, args=(data, w))
    print result[1]
    return result[0],result[0][0]+result[0][1]*np.sin(result[0][2]*data[0]+result[0][3]*np.pi)

def secondFit(data,result,tol=.1):
    def f(vars,data,weight,val):
        x = data[0]
        yaxis = data[1]       
        h = val[0]        
        a = val[1]
        f = val[2]
        p = vars[0]        
        model = (h+a*np.sin(f*x+p*np.pi))
        return (yaxis-model)*weight
    
    w = np.zeros(len(data[0]))
    i = 0
    for i in range(len(data[0])):
        d = data[1][i]
        y = result[0]+result[1]*np.sin(result[2]*data[0][i]+result[3]*np.pi)
        if np.abs(d-y)>tol:
            w[i] = 1
            
    val = result
    p = [result[3]-0.25]
    result = leastsq(f, p, args=(data, w,val))
    print result[1]
    return [val[0],val[1],val[2],result[0][0]],val[0]+val[1]*np.sin(val[2]*data[0]+result[0][0]*np.pi),w
 
def quartRes(filenum,xL,xH,side='r',tol=.1,n=100):
    data = getPrune(filenum)
    dat = getData(filenum)
    d,v = primeQuart(data,xL,xH)
    rf,yf = firstFit(d,v,n)
    rs,ys,w = secondFit(d,rf,tol)    
    lg = search(dat,xL)
    rg = search(dat,xH)
    vg = int(rg+lg)/2
    if side == 'r':
        dphi = (rs[3]-rf[3])%2
    else:
        dphi = (rf[3]-rs[3])%2
    return np.array([rf[2],rf[1],rf[3]*np.pi,rs[3]*np.pi,rf[0],0,1,dat[0][lg],dat[0][rg],dphi,vg,dat[0][vg]])
    