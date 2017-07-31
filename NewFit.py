# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 16:55:04 2014

@author: Simas
"""

import numpy as np
import pylab
import scipy
import random
from scipy.optimize import leastsq

#Function to turn the igor data files in to usable arrays
#0 is the gate voltage
#1 is the diagonal resistance
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

#Two search functions to retrieve index from arrays 
    
#function to be used when the gate voltage is being swept up to 0
#takes data and a value and generates closest index
def upsweep_search(seq, t):
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

#function to be used when the gate voltage is being swept down from 0
#takes data and a value and generates closest index
def downsweep_search(seq, t):
    min = 0
    max = len(seq) - 1
    while True:       
        m = (min + max) // 2        
        if max < min:
            if np.abs(t-seq[m]) < np.abs(t-seq[m+1]):          
                return m
            else:
                return m+1       
        if seq[m] > t:
            min = m + 1
        elif seq[m] < t:
            max = m - 1
        else:
            return m 

#Data visiualization and manipulation tools:

#Creates a small data segment from the initial data around Vp0 with a radius of dV  
def findValue(datas,Vp0,dV):
    min = 0
    max = 0    
    if datas[0,-1]<datas[0,0]:
        min = downsweep_search(datas[0],Vp0+dV)
        max = downsweep_search(datas[0],Vp0-dV)
    else:
        min = upsweep_search(datas[0],Vp0-dV)
        max = upsweep_search(datas[0],Vp0+dV)
    segment = np.ndarray(shape=(datas[0:2,min:max].shape))
    segment[0] = datas[0,min:max]
    segment[1] = datas[1,min:max]
    return segment
    
#Plots the data along with the fits for a given set of data and parameters
def plotVals(data,results,Vp0,dV):
    index = 0    
    if results[0,0]>results[-1,0]:
        index = downsweep_search(results[:,0],Vp0)
    else:
        index = upsweep_search(results[:,0],Vp0)
    real = findValue(data,results[index,0],dV)
    x = np.arange(len(real[1]))
    print results[index,7]    
    pylab.plot(real[0],real[1])
    pylab.plot(real[0],results[index,8]+results[index,1]+results[index,2]*x+results[index,3]*np.cos(results[index,4]*x-(results[index,5])*np.pi))
    pylab.plot(real[0],results[index,8]+results[index,1]+results[index,2]*x+results[index,3]*np.cos(results[index,4]*x-(results[index,6])*np.pi))
    pylab.show()
    
#Removes the average
def removeAverage(data):
    new = np.ndarray(shape=(data.shape))
    new[0] = data[0]
    new[1] = data[1] - sum(data[1])/len(data[1])
    return new

#Functions used in the least squares:

#Generates the initial guess based on the data
def generateVals(data):
    amplitude = np.std(data)
    frequency = 2*np.pi/1500
    p1 = random.gauss(0,2)
    p2 = random.gauss(0,2)
    height = sum(data)/len(data)
    slope = amplitude/len(data)
    a = np.arange(6,dtype=('float'))
    a[0] = height
    a[1] = slope    
    a[2] = amplitude
    a[3] = frequency
    a[4] = p1
    a[5] = p2    
    return a

#Step function with an exponential cut off, 100 units for 1/e           
def theta(x):
    w = 100
    result = np.arange(len(x),dtype=('float'))
    count = 0
    for i in x:
        if i >= 0:
            result[count]= 1
        else:
            result[count]= np.exp(i/w)
        count += 1
    return np.asarray(result)
    
#Creates a weight of either one or zero which is used to simply the least sq fit        
def weight(data,points,count):
    new = np.ndarray(shape=(data.shape))
    new[1] = 1
    new[0] = data[0]
    j = points   
    while j < (len(new[1])-points):
        counts = 0        
        for k in range(1,points+1):
            if np.abs(data[1,j]-data[1,j+k]) > .1:
                counts += 1
            if np.abs(data[1,j]-data[1,j-k]) > .1:
                counts += 1
        if counts > count:
            new[1,j] = 0
        j += 1
    return new    
    
#The least squares funtion with a weight
def testLS(vars,DATA,eps,weight):
    xaxis = DATA[0]
    yaxis = DATA[1]    
    Vp0 = len(xaxis)/2     
    a = vars[2]
    f = vars[3]
    pp = vars[4]
    pm = vars[5]
    h = vars[0]
    l = vars[1]
    x = np.arange(len(yaxis))
    model = theta(Vp0-x)*(h+l*x+a*np.cos(f*x-pp*np.pi))
    model += theta(x-Vp0)*(h+l*x+a*np.cos(f*x-pm*np.pi))
    return (yaxis-model)*weight[1]/eps

#The main function to run, fits slips to the entire data file
#Saves the results into a numpy file
def mainRun(name,points,dV):
    data=parse(name)        #create the data array
    weights = weight(data,50,25)  #generate leastsq weight
    min = 0 #generate lower index for search, to be incremented for total searches
    max = 0 #generate upper index for search
    if data[0,-1]<data[0,0]: #assign values to the min and max index
        min = downsweep_search(data[0],data[0,0]-dV)
        max = downsweep_search(data[0],data[0,-1]+dV)
    else:
        min = upsweep_search(data[0],data[0,0]+dV)
        max = upsweep_search(data[0],data[0,-1]-dV)
    step = int((max-min)/(points)) #creates step length for index based on required searches
    results = np.ndarray(shape=(points+1,9)) #creates a result array with the proper shape
    i = 0
    print 'Starting loops'
    while i <= points:       
        results[i,0] = data[0,min] #stores Vp0 value      
        segment = findValue(data,data[0,min],dV) #creates data segment to be used in the fit
        weightsegment = findValue(weights,weights[0,min],dV)#creates weight segment to be used
        averageheight = sum(segment[1])/len(segment[1]) #centers the resistance around 0 for proper fitting and weighing       
        segment[1] -= averageheight   #remove average for theta function     
        result = leastsq(testLS, generateVals(segment[1]), args=(segment, 1,weightsegment))
        if result[1]<5:       
            phi1 = result[0][4]
            phi2 = result[0][5]
            ans = phi1 - phi2        
            while (ans)<0:
                ans+=2
            while(ans)>2:
                ans-=2        
            results[i,1] = result[0][0]
            results[i,2] = result[0][1]
            results[i,3] = result[0][2]
            results[i,4] = result[0][3]
            results[i,5] = phi1
            results[i,6] = phi2
            results[i,7] = ans
            results[i,8] = averageheight
        min += step
        i += 1
    np.save(name+str(dV),results)
    return results

#The commands to run
trytest = mainRun('rd20a.vp',120,7)
pylab.scatter(trytest[:,0],trytest[:,7])
pylab.show()


    