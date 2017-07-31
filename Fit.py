import numpy as np
import pylab
import scipy.io
import scipy.optimize
import random
from scipy.optimize import leastsq

#Parses the data file into an array
#0 is the gate voltage
#1 is the diagonal resistance
def parse(filenum):
    f = open('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\rd'+str(filenum)+'a.vp', 'r')
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

#Removes a slope from the data for better graphing
def cleanData(filenum):
    data = parse(filenum)
    def minFunc(a):
        return max(data[1] + a*data[0]) - min(data[1] + a*data[0])
    out = scipy.optimize.fmin(minFunc,0,disp=False)
    data[1]=data[1]+out[0]*data[0]
    return data

#Removes the mean so the values are centered around 0
def removeAverage(data):
    new = np.ndarray(shape=(data.shape))
    new[0] = data[0]
    new[1] = data[1] - sum(data[1])/len(data[1])
    return new,sum(data[1])/len(data[1])

#Creates an array of 1s for all good data and 0 for bad data
#takes width of points and count for good data as args, default 50 and 25 resp
#Returns array with same x axis, y values are 1s and 0s
def prune(data,points=50,count=25):
    new = np.ndarray(shape=(data.shape))
    new[1] = 1
    new[0] = data[0]
    j = points   
    while j < (len(new[1])-points):
        counts = 0        
        for k in range(1,points+1):
            if np.abs(data[1,j]-data[1,j+k]) > .036:
                counts += 1
            if np.abs(data[1,j]-data[1,j-k]) > .036:
                counts += 1
        if counts > count:
            new[1,j] = 0
        j += 1
    return new
#Generates a smooth step function with exp decay of width w
def theta(x):
    w = 1
    out = list(x)
    for i in out:
        if i >= 0:
            i =1
        else:
            i = np.exp(i/w)            
    return np.asarray(out)
#Least squares function
def LS(var,data,Vp0,eps,prune):
    xaxis = data[0]
    yaxis = data[1]    
    Vp0 = Vp0   
    a = var[1]
    f = var[0]
    pl = var[2]
    pr = var[3]
    h = var[4]
    l = var[5]
    model = theta(Vp0-xaxis)*(h+l*xaxis+a*np.cos(f*xaxis-pr*np.pi))
    model += theta(xaxis-Vp0)*(h+l*xaxis+a*np.cos(f*xaxis-pl*np.pi))
    return (yaxis-model)*prune[1]/eps
#Initial Guess
def generateVals(data):
    amplitude = np.std(data[1])
    frequency = 2*np.pi/10
    pl = random.gauss(0,2)
    pr = random.gauss(0,2)
    height = sum(data[1])/len(data[1])
    slope = amplitude/len(data[1])
    a = np.arange(6,dtype=('float'))
    a[4] = height
    a[5] = slope    
    a[1] = amplitude
    a[0] = frequency
    a[2] = pl
    a[3] = pr    
    return a

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
            
def findStep(data,steps):
    minv = data[0,0] + 0.05
    maxv = data[0,-1] - 0.05
    mini = search(data,minv)
    maxi = search(data,maxv)
    step = int((maxi-mini)/steps)
    return step
    
def bounds(data,Vp0,dV=7):
    minv = data[0,0] + 0.05
    maxv = data[0,-1] - 0.05
    mini = search(data,minv)
    maxi = search(data,maxv)
    lbound = 0
    rbound = 0
    if (Vp0-dV)<minv:
        lbound = mini
    else:
        lbound = search(data,(Vp0-dV))
    if (Vp0+dV)>maxv:
        rbound = maxi
    else:
        rbound = search(data,(Vp0+dV))
    return lbound,rbound
    
def main(filenum,steps=120,dV=7 ):
    data = cleanData(filenum)
    out = np.ndarray(shape=(steps,12))
    pruned = prune(data)
    nomean,out[:,4] = removeAverage(data)
    step = findStep(data,steps)
    li = 0
    ri = 0
    vpi= 0
    vp = 0
    i = 0
    while i<steps:
        vpi = step+i*step
        vp = nomean[0][vpi]
        li,ri = bounds(nomean,vp)
        out[i,7] = li
        out[i,8] = ri
        out[i,10] = vpi
        out[i,11] = vp
        tempd = nomean[:,li:ri+1]
        tempp = pruned[:,li:ri+1]
        result = leastsq(LS, generateVals(tempd), args=(tempd,vp,.01,tempp))
        if result[1]<5:
            out[i,0] = result[0][0]
            out[i,1] = result[0][1]
            out[i,2] = result[0][2]
            out[i,3] = result[0][3]
            out[i,4] += result[0][4]
            out[i,5] = result[0][5]
            out[i,6] = quality(tempd,out[i])
            out[i,9] = (result[0][2]-result[0][3])%2
        i += 1
    return out

#Returns Lorentzian of fit quality, maybe use pruned data for fit         
def quality(data,var):
    xaxis = data[0]
    yaxis = data[1]    
    Vp0 = var[11]   
    a = var[1]
    f = var[0]
    pl = var[2]
    pr = var[3]
    h = var[4]
    l = var[5]
    yr = (h+l*xaxis+a*np.cos(f*xaxis-pr*np.pi))
    yl = (h+l*xaxis+a*np.cos(f*xaxis-pl*np.pi))
    errl = theta(xaxis-Vp0)/(1+((yl-yaxis)/.01)**2)
    errr = theta(Vp0-xaxis)/(1+((yr-yaxis)/.01)**2)
    qual = sum(errl)+sum(errr)
    return qual
    
    
    