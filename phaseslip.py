# -*- coding: utf-8 -*-
"""
Created on Mon Apr 07 12:38:49 2014

@author: Simas
"""
import numpy as np
import pylab
import scipy
from scipy.optimize import curve_fit

def parse(filename):
    f = open('C:\\Users\\Simas\\Desktop\\Research\\Phase Slip Analysis\\phase slip data\\'+filename, 'r')
    raw = f.readline()
    parse = raw.split('\r',-1)
    data =np.ndarray(shape=(2,len(parse)-2))   
    count = 0
    for a in parse: 
        if count > 0 and count < (len(parse)-1):                               
            data[0,count-1]=float(a.split('\t',-1)[0])
            data[1,count-1]=float(a.split('\t',-1)[1])
        count += 1    
    return data
    
def removeSlope(data):
    line = np.polyfit(data[0],data[1],1)
    linefn = np.poly1d(line)
    new =np.ndarray(shape=(data.shape))
    new[0] = data[0]
    new[1] = data[1]-linefn(data[0])
    return new
        
    
def smoothData(data,counts):
    n=counts
    new =np.ndarray(shape=(data.shape))
    new[0] = data[0]
    data = data[1]
    f = len(data)
    clean = np.zeros(f)
    for i in range(n,f-n):
        total = data[i]
        count = 1
        for j in range(1,n):
            if (data[i+j] - data[i]) < .1 and  (data[i+j] - data[i]) > -.1:          
                total += data[i+j]
                count +=1
            if (data[i-j] - data[i]) < .1 and  (data[i-j] - data[i]) > -.1:           
                total += data[i-j]
                count +=1    
            clean[i] = total/count
    new[1] = clean
    return new

def removeSpike(data,points):
    new = np.ndarray(shape=(data.shape))
    new[0]=data[0]
    new[1]=data[1]
    j = points
    while j>0:
        for i in range(1,len(data[1])-j):
            if ((new[1,i]-new[1,i-j]) >.02 or (new[1,i]-new[1,i-j]) <-.02) and ((new[1,i]-new[1,i+j]) >.02 or (new[1,i]-new[1,i+j]) <-.02):
                new[1,i] = (new[1,i-j]+new[1,i+j])/2
        j -= 1
    return new
def goodData(name,smoothcount,spikecount):
    data = parse(name)
    slope = removeSlope(data)
    smooth= smoothData(slope,smoothcount)
    nospike = removeSpike(smooth,spikecount)
    return nospike
#test function for wave
def func(x,a,b,c,d,f):
    return a+b*x+c*np.cos(d*x-f)
    
def func(x,a,b,c,d):
    return a+b*np.cos(c*x-d)
data = parse('rd20a.vp')
slope = removeSlope(data)
smooth= smoothData(slope,50)
nospike = removeSpike(smooth,5)
nospike2 = removeSpike(nospike,4)
nospike3 = removeSpike(nospike2,3)
nospike4 = removeSpike(nospike3,2)
nospike5 = removeSpike(nospike4,1)
pylab.subplot(5,1,1)
pylab.plot(nospike[0],nospike[1])
pylab.subplot(5,1,2)
pylab.plot(nospike2[0],nospike2[1])
pylab.subplot(5,1,3)
pylab.plot(nospike3[0],nospike3[1])
pylab.subplot(5,1,4)
pylab.plot(nospike4[0],nospike4[1])
pylab.subplot(5,1,5)
pylab.plot(nospike5[0],nospike5[1])
pylab.show()

popt, pcov = curve_fit(func, test[0][6775:7770],test[1][6775:7770])

transform = np.abs(np.fft.rfft(test[1][4213:7870]))**2
field = (test[0][4213]-test[0][7870])*1.
xaxis = np.arange(len(transform))/(field)
transform[0]=0
pylab.plot(xaxis,transform)

slope = removeSlope(data)
pylab.plot(slope[0][4200:8040],slope[1][4200:8040])

pylab.plot(test[0][4200:8040],func(test[0][4200:8040],test2[0],test2[1],test2[2],-.43))

pylab.plot(test[0][4200:8040],func(test[0][4200:8040],test2[0],test2[1],test2[2],test2[3]))
test2, pcov = curve_fit(func,test[0][4800:6450],test[1][4800:6450])