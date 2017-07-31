import matplotlib
import numpy as np
import pylab
import scipy.io
import scipy.optimize
import Histogram as H
#load the matlab fit file
#mat = scipy.io.loadmat('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\output.mat')
#results = mat['result']
#Parses the data file into an array
#0 is the gate voltage
#1 is the diagonal resistance
def parseFit(filenum):
    f = open('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\'+str(filenum), 'r')
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
    return data


def parse(filenum):
    f = open('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\'+str(filenum), 'r')
    raw = f.readline()
    parse = raw.split('\r',-1)
    data =np.ndarray(shape=(2,len(parse)-2))   
    count = 0
    for a in parse: 
        if count > 0 and count < (len(parse)-1):                               
            data[0,count-1]=float(a.split('\t',-1)[0])
            data[1,count-1]=float(a.split('\t',-1)[1])
        count += 1    
    if data[0,0]>data[0,-1]:
        data = np.fliplr(data)    
    return data

def getResults(unpruned=True,n=0):
    mat = scipy.io.loadmat('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\output.mat')
    results = mat['result']    
    if unpruned == True:
        return results
    else:
        index = 0        
        for i in results:
            data = parse(n*10+index)
            pruned = parsePrune(n*10+index)             
            for j in i:
                j[7] = list(data[0]).index(pruned[0][j[7]-1])
                j[8] = list(data[0]).index(pruned[0][j[8]-1])
                j[10] = list(data[0]).index(pruned[0][j[10]-1])
            index += 1
        return results
def parsePrune(filenum):
    f = open('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\'+str(filenum)+'x', 'r')
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
    return data
#give regular file number name and calculate the width in vstep of slips
def slipWidth(results,filenum,n):
    slips = results[filenum-n*10,:,9]
    out = list()
    start = 0
    stop = 0
    check = 0
    i= 0
    while(i<119):
        if (check == 0 and (slips[i]>0.20 and slips[i]<1.80)):          
            start = i
            check = 1
        if np.abs(slips[i]-slips[i+1])<0.05:
            i = i + 1
        elif check>0 and np.abs(slips[i]-slips[i+1])>0.05:
            stop = i
            if start != stop:
                out.append([start,stop])
            check = 0
            i = i +1
        else:
            i = i + 1
    return out
#finds slip with greatest fit and returns a list with index, slip value, and fit quality           
def findMaxFit(results,filenum,widths,n):
    filenum = filenum - n*10    
    index = 0
    final = list()
    fit = results[filenum,:,6]
    for out in widths:
        big = 0
        i = out[0]
        while i <= out[1]:
            if fit[i] > big:
                big = fit[i]
                index = i
            i = i + 1        
        final.append([index,results[filenum,index,9],results[filenum,index,6]])
    return final

#Constructs a list filled with MaxFits for all sweeps
def findAllMaxFit(results,n):
    out = list()
    for i in range(len(results)):
        width = slipWidth(results,i+n*10,n)
        out.append(findMaxFit(results,i+n*10,width,n))
    return out

#returns threshold for quality, 1.5 std above mean
def returnThresh(results,x,filenum=-1,n=0):
    out = list()
    maxFit = list()
    maxFit = findAllMaxFit(results,n)
    if filenum==-1:
        for i in maxFit:
            for j in i:
                out.append(j[2])
    else:
        for j in maxFit[filenum-n*10]:
            out.append(j[2])
    return np.mean(out)+x*np.std(out)
    
#returns an array of fits that meet a threshold value and are outisde a 
#redundancy limit, returns fit step, delta phi, and quality    
def getBestFit(MaxFit,threshold):   
    fits = list(MaxFit)     
    out = list()
    for i in fits:
        if i[2]>threshold:
            out.append(i)
    hold = list()
    out = list(out)
    j = 0
    index = 0
    final = list() 
    while j <12:
        if index < len(out) and (out[index][0] >= j*10 and out[index][0] < (j*10+9)):
            hold.append(out[index])
            index += 1
        elif len(hold)==0:
            j +=1
        else:
            hold = np.array(hold)
            best = np.max(hold[:,2])
            x = list(hold)[list(hold[:,2]).index(best)]
            final.append(list(x))
            hold = list()
            j +=1
    return final

def getBestFit2(MaxFit,threshold):
    fits = list(MaxFit)     
    out = list()
    for i in fits:
        if i[2]>threshold:
            out.append(i)
    step = 0
    hold = list()
    final = list()
    while len(out) > 0:
        if out[0][0] < step + 3:
            hold.append(out[0])
            out.remove(out[0])
            if len(out)==0:
                hold = np.array(hold)
                best = np.max(hold[:,2])
                x = list(hold)[list(hold[:,2]).index(best)]
                final.append(list(x))
                hold = list()
        else:
            if len(hold) > 0:
                hold = np.array(hold)
                best = np.max(hold[:,2])
                x = list(hold)[list(hold[:,2]).index(best)]
                final.append(list(x))
                hold = list()
            step = out[0][0]
    return final
    
#Returns best fits for selected data
def fit(results,filenum,n,threshold):
    width = slipWidth(results,filenum,n) #find the nontrivial fit widths
    mf = findMaxFit(results,filenum,width,n) #find the fits with best quality
    bf = getBestFit2(mf,threshold)
    return bf
#Plots the left and right fit 
def plotfit(results,filenum,n,step,back=True,max=False):
    xa = results[filenum-n*10,step,:]    
    freq = xa[0]
    amp = xa[1] 
    phasel = xa[2]
    phaser = xa[3]
    height = xa[4]
    slope = xa[5]
    posleft = xa[7] - 1
    posright = xa[8] - 1
    pos =  xa[10] -1
    v0 = xa[11]
    xposl = list()
    xposr = list()
    data = parse(filenum)        
    xposl = data[0,posleft:pos+1]
    xposr = data[0,pos:posright+1]
    if max == True:
        xposl = data[0,0:pos+1]
        xposr = data[0,pos:]
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

#plot a line with a given phase difference
def plotline(results,filenum,step,n,phi):
    xa = results[filenum-n*10,step,:]    
    freq = xa[0]
    amp = xa[1] 
    phasel = xa[2]
    height = xa[4]
    slope = xa[5]
    xpos = list()
    data = parse(filenum)        
    xpos = data[0]
    pos =  xa[10] -1      
    ytarg = height + amp*np.sin(freq*xpos + phasel + phi) + slope*(xpos - data[0,pos])
    pylab.plot(xpos,ytarg,'go',markersize=4,mew=0)


#Removes a slope from the data for better graphing
def cleanData(filenum):
    data = parse(filenum)
    def minFunc(a):
        return max(data[1] + a*data[0]) - min(data[1] + a*data[0])
    out = scipy.optimize.fmin(minFunc,0,disp=False)
    data[1]=data[1]+out[0]*data[0]
    return data
#Given an x value from the data, returns a list of y values for fits within
#the bounds
#Return in a form (fit number,y value,phase value) for left and then right    
def whichSlip(x,data,results,filenum,bestFits,n):       
    bfindex = np.array(bestFits)[:,0]
    filenum = filenum - n*10    
    hold = list()
    out = list()
    x = data[0,x]
    for i in bfindex:
        if x >= results[filenum,i,7] and x <= results[filenum,i,8]:
            hold.append(i)
    for i in hold:
        xa = results[filenum,i,:]    
        freq = xa[0]
        amp = xa[1] 
        phasel = xa[2]
        phaser = xa[3]
        height = xa[4]
        slope = xa[5]
        v0 =  xa[11]
        yl = height + amp*np.sin(freq*x + phasel) + slope*(x - v0)
        yr = height + amp*np.sin(freq*x + phaser)  + slope*(x - v0)        
        out.append([i,yl,phasel,xa[9]])
        out.append([i,yr,phaser,xa[9]])
    return out
    
def whichSlip2(x,data,results,filenum,bestFits,n):       
    bfindex = np.array(bestFits)[:,0]
    filenum = filenum - n*10    
    hold = list()
    out = list()
    for i in bfindex:
        if x >= results[filenum,i,7]-1 and x <= results[filenum,i,8]-1:
            hold.append(i)
    for i in hold:
        xa = results[filenum,i,:]    
        freq = xa[0]
        amp = xa[1] 
        phasel = xa[2]
        phaser = xa[3]
        height = xa[4]
        slope = xa[5]
        pos =  xa[10] -1
        xpos = data[0,x]
        yl = height + amp*np.sin(freq*xpos + phasel) + slope*(xpos - data[0,pos])
        yr = height + amp*np.sin(freq*xpos + phaser)  + slope*(xpos - data[0,pos])
        if np.abs(yl-yr)>0.02:       
            out.append([i,yl,phasel])
            out.append([i,yr,phaser])

    return out
#Creates an array of 1s for all good data and 0 for bad data
#takes width of points and count for good data as args, default 50 and 25 resp
#Returns array with same x axis, y values are 1s and 0s
def prune(data,points=30,count=20,tol=0.072):
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


#The meat of the analysis, goes through the data point by point and checks if the
#point is real or noise. If it's real consults whichSlip and checks to see if 
#there is an appropriate fit for the point. If a fit exists checks adjacent point
#for a fit and phase change. If there is a phase change to an appropriate fit 
#records the value.
def slipCount(result,filenum,n,threshold=-1,points=30,count=20,tol=0.072,minv=0,maxv=0):
    datas = cleanData(filenum) #get the data
    results = np.array(result)
    width = slipWidth(results,filenum,n) #find the nontrivial fit widths
    mf = findMaxFit(results,filenum,width,n) #find the fits with best quality
    if threshold == -1:
        threshold = returnThresh(results,0,filenum,n) #Mean quality of fits becomes threshold
    bf = getBestFit2(mf,threshold) #clean the fits for some threshold
    data = np.array(prunedData(datas,points,count,tol)) #Try using pruned data instead of checking points
    final = list()
    hold0 = list()
    hold1 = list()
    stepval = list()   
    my = 0.02
    results[filenum-n*10] = resValue(results,datas,filenum,n)#takes real data and turns indeces into values
    if minv == 0:
        minv = data[0,0] + 0.05
    if maxv == 0:
        maxv = data[0,-1] - 0.05
    for i in range(len(data[0])-1):#runs through the data
        if data[0,i] > minv and data[0,i] <maxv:
            pass
        else:
            continue
        if len(bf)>0:   #Checks to see if there are any fits for this data set         
            hold0 = list(whichSlip(i,data,results,filenum,bf,n)) #initializes the fits to check
        else:
            break
        best = list()
        miny = 1
        left = 0
        if len(hold0)!=0: #check for fits at this point
            miny = 1 #set initial difference
            index = 0
            bindex = 0
            for j in hold0: #cycle through fits at this point to minimize diff
                if np.abs(data[1,i]-j[1]) < miny:
                    bindex = index
                    miny = np.abs(data[1,i]-j[1])
                    left = index%2
                index += 1
            best = list(hold0[bindex]) #store the info for this fit
        if miny < my: #noise check for first point
            hold1 = list(whichSlip(i+1,data,results,filenum,bf,n)) #get fits for next point
            if len(hold1)!=0:
                miny = 1
                index = 0
                bindex = -1
                for j in hold1:
                    if best[0]==j[0] and best[1]!=j[1] and np.abs(data[1,i+1]-j[1]) < miny: #forces the code to only check compatible fits 
                        bindex = int(index)
                        miny = np.abs(data[1,i+1]-j[1])
                    index += 1
                if miny < my and best[2]!=hold1[bindex][2]: #second noise check                  
                    if left ==0:
                        fin = ((best[2]-hold1[bindex][2])/np.pi)%2.0
                    else:
                        fin = ((hold1[bindex][2]-best[2])/np.pi)%2.0 #clean up the slip 
                    final.append(fin)
                    if best[0] not in stepval:
                        print best[0],best[3]
                        stepval.append(best[0])
    return final


#same as slipCount with a min an max voltage window
#if no min or max, +1 mV and -1 mV from ends    
def slipCountMM(results,filenum,threshold,n,minv=0,maxv=0):
    data = cleanData(filenum) 
    width = slipWidth(results,filenum,n)
    mf = findMaxFit(results,filenum,width,n)
    bf = getBestFit2(mf,threshold)
    final = list()
    hold0 = list()
    hold1 = list()
    my = 0.02
    if minv ==0 and maxv==0:
        minv = data[0,0] + 0.05
        maxv = data[0,-1] - 0.05
    for i in range(len(data[0])-1):
        if data[0,i] > minv and data[0,i] <maxv:
            pass
        else:
            continue
        if len(bf)>0:            
            hold0 = list(whichSlip(i,data,results,filenum,bf,n))
        best = list()
        miny = 1
        left = 0
        if len(hold0)!=0:
            miny = 1
            index = 0
            bindex = 0
            for j in hold0:
                if np.abs(data[1,i]-j[1]) < miny:
                    bindex = index
                    miny = np.abs(data[1,i]-j[1])
                    left = index%2
                index += 1
            best = list(hold0[bindex])
        if miny < my:
            hold1 = list(whichSlip(i+1,data,results,filenum,bf,n))
            if len(hold1)!=0:
                miny = 1
                index = 0
                bindex = -1
                for j in hold1:
                    if best[0]==j[0] and best[1]!=j[1] and np.abs(data[1,i+1]-j[1]) < miny:
                        bindex = int(index)
                        miny = np.abs(data[1,i+1]-j[1])
                    index += 1
                if miny < my and best[2]!=hold1[bindex][2]:                   
                    if left ==0:
                        fin = ((best[2]-hold1[bindex][2])/np.pi)%2.0
                    else:
                        fin = ((hold1[bindex][2]-best[2])/np.pi)%2.0                    
                    final.append(fin)
    return final
#Uses pruned data, main run!    
def slipCountMM2(results,filenum,threshold,n,minv=0,maxv=0):
    data = cleanData(filenum) 
    width = slipWidth(results,filenum,n)
    mf = findMaxFit(results,filenum,width,n)
    bf = getBestFit2(mf,threshold)
    pruned = prune(data) 
    final = list()
    hold0 = list()
    hold1 = list()
    my = 0.02
    if minv ==0 and maxv==0:
        minv = data[0,0] + 0.05
        maxv = data[0,-1] - 0.05
    for i in range(len(data[0])-1):
        if data[0,i] > minv and data[0,i] <maxv:
            pass
        else:
            continue
        if pruned[1,i]==1 and pruned[1,i+1]==1:
            pass
        else:
            continue
        if len(bf)>0:            
            hold0 = list(whichSlip2(i,data,results,filenum,bf,n))
        best = list()
        miny = 1
        left = 0
        if len(hold0)!=0:
            miny = 1
            index = 0
            bindex = 0
            for j in hold0:
                if np.abs(data[1,i]-j[1]) < miny:
                    bindex = index
                    miny = np.abs(data[1,i]-j[1])
                    left = index%2
                index += 1
            best = list(hold0[bindex])
        if miny < my:
            hold1 = list(whichSlip2(i+1,data,results,filenum,bf,n))
            if len(hold1)!=0:
                miny = 1
                index = 0
                bindex = -1
                for j in hold1:
                    if best[0]==j[0] and np.abs(data[1,i+1]-j[1]) < miny:
                        bindex = int(index)
                        miny = np.abs(data[1,i+1]-j[1])
                    index += 1
                if miny < my and best[2]!=hold1[bindex][2]:                   
                    if left ==0:
                        fin = ((best[2]-hold1[bindex][2])/np.pi)%2.0
                    else:
                        fin = ((hold1[bindex][2]-best[2])/np.pi)%2.0                    
                    final.append(fin)
    return final
    
#Runs slipCount for all the files from a results array
#Uses pruned data and whichSlip 
def total(results,n,threshold=-1,points=30,count=20,tol=0.072,minv=0,maxv=0):
    final = list()
    out = list()
    for i in range(len(results)):
        final.append(slipCount(results,i+n*10,n,threshold,points,count,tol,minv,maxv))
    for i in final:
        for j in i:
            out.append(j)
    return out
    
    
#same as total with a min and max voltage window
def totalMM(results,threshold,n,minv=0,maxv=0):
    final = list()
    out = list()
    for i in range(len(results)):
        final.append(slipCountMM(results,i+n*10,threshold,n,minv,maxv))
    for i in final:
        for j in i:
            out.append(j)
    return out
    
#Full fit, uses prune values and can set edge which are default at 0    
def totalMM2(results,threshold,n,minv=0,maxv=0):
    final = list()
    out = list()
    for i in range(len(results)):
        final.append(slipCountMM2(results,i+n*10,threshold,n,minv,maxv))
    for i in final:
        for j in i:
            out.append(j)
    return out

#Create two column matrix for exporting
def twoCol(data):
    out = [[0] * 2 for i in range(len(data[0]))]
    for i in range(len(data[0])):
        out[i][0] = data[0][i]
        out[i][1] = data[1][i]
    return out
    
def prunedData(data,points=30,count=20,tol=0.072):
    pruned = prune(data,points,count,tol)
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

def returnPruned(filename,points=30,count=20,tol=0.072):
    data = parse(filename)
    pruned = prunedData(data,points,count,tol)
    out = twoCol(pruned)
    np.savetxt(str(filename)+'x',out)

#Plots histogram using output of fit
def plot(out):
    bins = np.arange(0,2.0,.01)
    matplotlib.pyplot.hist(out,bins)     
#bins = np.arange(0,2.0,.01)
#matplotlib.pyplot.hist(out,bins)
def fakeTele(y1,y2,start,stop):
    hold = np.array(y1)    
    line = 0 
    i = start
    while i <= stop:
        width = int(np.random.rand()*500+1)
        if line ==1:
            hold[i:i+width+1]=np.array(y2[i:i+width+1])
            line = 0
        else:
            line = 1
        i = i + width
    return hold

def resValue(results,data,filenum,n):
    i = results[filenum-n*10]
    for j in i:
        j[7] = data[0][j[7]-1]
        j[8] = data[0][j[8]-1]
    return i

fs = findSlips(filenum,dV,pruned,minv,maxv,mindiff)
alll = list()
for i in fs:
    for j in i:
        alll.append(j)