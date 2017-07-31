import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize


#Initialzie local data storage
data = list()
prune = list()
results = list()
bestFits = list()


def init():
    parseData()
    parsePrune('soft')
    getPrune()
    return getData(), getPrune(), [12.145, 0.09, 2 * np.pi / 7.0]


# set and get functions for data, prune, and results
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


# initialzie data,prune,and results
def parseData():
    f = open('filenames')
    raw = f.read()
    parse = raw.split('\n', -1)
    out = list()
    f.close()
    for i in parse:
        if i == '':
            continue
        f = open('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\' + str(i), 'r')
        raw = f.readline()
        parse = raw.split('\r', -1)
        data = np.ndarray(shape=(2, len(parse) - 2))
        count = 0
        for a in parse:
            if count > 0 and count < (len(parse) - 1):
                hold = a.split('\t', -1)
                data[0, count - 1] = float(hold[0])
                data[1, count - 1] = float(hold[1])
            count += 1
        if data[0, 0] > data[0, -1]:
            data = np.fliplr(data)

        # clean the data
        def minFunc(a):
            return max(data[1] + a * data[0]) - min(data[1] + a * data[0])

        slope = scipy.optimize.fmin(minFunc, 0, disp=False)
        data[1] = data[1] + slope[0] * data[0]
        out.append(data)
    setData(out)


def parsePrune(filetype):
    f = open('prunenames')
    raw = f.read()
    parse = raw.split('\n', -1)
    out = list()
    f.close()
    for i in parse:
        f = open('C:\\Users\\Simas\Desktop\\Simon Matlab\\matlab\\' + filetype + '\\' + str(i), 'r')
        raw = f.read()
        parse = raw.split('\n', -1)
        data = np.ndarray(shape=(2, len(parse) - 1))
        count = 0
        for a in parse:
            if count < (len(parse) - 1):
                data[0, count] = float(a.split(' ', -1)[0])
                data[1, count] = float(a.split(' ', -1)[1])
            count += 1
        if data[0, 0] > data[0, -1]:
            data = np.fliplr(data)
        out.append(data)
    setPrune(out)


def bruteforce(data, params):
    x = np.array(data[0])
    y = np.array(data[1])
    best = list()
    cost = np.inf
    left = step(data, 'left')
    right = step(data, 'right')
    h = params[0]
    a = params[1]
    f = params[2]
    H = np.linspace(h-h*0.05,h+h*0.05,20)
    A = np.linspace(a - a * 0.1, a + a * 0.1, 20)
    F = np.linspace(f - f * 0.1, f + f * 0.1, 20)
    PL = np.linspace(0, 2, 20)
    PR = np.linspace(0, 2, 20)
    for hi in H:
        for fi in F:
            for ai in A:
                for pli in PL:
                    for pri in PR:
                        c = np.sum(((y - (hi + ai * np.sin(fi * x + pli * np.pi)))*left) ** 2
                                   + ((y - (hi + ai * np.sin(fi * x + pri * np.pi)))*right) ** 2)
                        if c < cost:
                            cost = c
                            best = [hi, ai, fi, pli, pri]
    return best


def step(data, direction):
    l = len(data[0])
    out = list()
    mid = data[0][l/2]
    if direction == 'left':
        for x in data[0]:
            if x > mid:
                if np.exp((mid-x)/.5) > 0.10:
                    out.append(np.exp((mid-x)/.5))
                else:
                    out.append(0.0)
            else:
                out.append(1.0)
    else:
        for x in data[0]:
            if x < mid:
                if np.exp((x-mid)/.5) > 0.10:
                    out.append(np.exp((x-mid)/.5))
                else:
                    out.append(0.0)
            else:
                out.append(1.0)
    return np.array(out)


def plotFig(data, fitbit, left=0, right=0):
    xa = fitbit
    freq = xa[2]
    amp = xa[1]
    phasel = xa[3]*np.pi
    phaser = xa[4]*np.pi
    height =  xa[0]
    if left == 0:
        posleft = data[0][0]
    else:
        posleft = left
    if right == 0:
        posright = data[0][-1]
    else:
        posright = right
    xpos = np.linspace(posleft,posright,10000)
    ytargl = height + amp*np.sin(freq*xpos + phasel)
    ytargr = height + amp*np.sin(freq*xpos + phaser)
    plt.plot(data[0],data[1])
    plt.plot(xpos,ytargl,'k',linewidth=2)
    if phaser != 0:
        plt.plot(xpos,ytargr,'r',linewidth=2)


def cost(data, params):
    h = params[0]
    a = params[1]
    f = params[2]
    pl = params[3]
    pr = params[4]
    x = np.array(data[0])
    y = np.array(data[1])
    left = step(data, 'left')
    right = step(data, 'right')
    c = np.sum(((y - (h + a * np.sin(f * x + pl * np.pi))) * left) ** 2
               + ((y - (h + a * np.sin(f * x + pr * np.pi))) * right) ** 2)
    return c

def sectiondata(data,dV=7.0,n=20):
    zm = np.array(data)
    xmin = zm[0][0]
    xmax = zm[0][-1]
    step = (xmax-xmin)/n
    cstep = xmin
    out = list()
    while cstep<xmax:
        if cstep-dV/2.0 < xmin:
            xl = xmin
        else:
            xl = cstep-dV/2.0
        if cstep+dV/2.0 > xmax:
            xr = xmax
        else:
            xr = cstep+dV/2.0
        li = search(zm,xl)
        ri = search(zm,xr)
        out.append([xl,xr,np.array(zm[:,li:ri])])
        cstep += step
    return out



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


def multirun(data,n):
    out = list()
    sd = sectiondata(data, n)
    for i in sd:
        d = i[2]
        left = step(d, 'left')
        right = step(d, 'right')
        v = val(d)
        out.append(bruteforce(d,v))
    return out

def val(data):
    h = np.mean(data[1])
    a = np.std(data[1])
    f = 2*np.pi/7.0
    return [h,a,f]