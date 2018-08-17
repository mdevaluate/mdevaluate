import numpy as np
import scipy
from mdevaluate.pbc import pbc_diff
import sys
import collections
import multiprocessing

def Hbondlist(don,H_,akz,maxangle=30,maxdistance=0.35):
    """
    Returns the hydrogen bonds of each donor-group for a frame in the trajectory.
    
    Args:
        don: The coordintes of the donor atoms of each donor-group.
        H_: The coordinates of the hydrogen atoms of each donor-group.
        akz: The coordinates of each acceptor atom.
        maxangle: The maximal angle between donor-group and acceptor
        maxdistance: The maximal distance between donor-atom and acceptor
    """
    
    angle = np.cos(maxangle*2*np.pi/360)
    don = don%don.box.diagonal()
    akz = akz%akz.box.diagonal()
    par = dict(balanced_tree = False, compact_nodes = False , leafsize = 32)
    tree = scipy.spatial.cKDTree(akz,**par,boxsize=akz.box.diagonal())
    na = tree.query(don,k = 10)

    a = pbc_diff(akz[na[1]],don[:,np.newaxis,:], box = akz.box.diagonal())
    b = pbc_diff(H_,don, box = akz.box.diagonal())

    scalar = (np.multiply(a,b[:,np.newaxis,:]).sum(2))
    mask = (na[0]>0)
    c = ((b**2).sum(1))**0.5
    fin = np.divide(scalar,na[0], where=mask)/c[:,np.newaxis]

    finalmask = (na[0]<maxdistance)& ((fin) > angle) & (na[0]>0)


    result = np.array([na[1][i, x] for i, x in enumerate(finalmask)])
    return result

def Veranderungsliste(donliste,Hliste,akzliste,donqty=1):
    """
    Returns the list of all hydrogenbond changes as a string of all changes of each donor-atom and the start configuration.
    
    Args:
        donliste: The subset of a Trajectory with the donor-atom of each donor-group.
        Hliste: The subset of a Trajectory with the hydrogen atom of each donor-group.
        akzliste: The subset of a Trajectory with all acceptor atoms.
    """
    mask=[]
    for i in range(len(donliste[0])):
        j = 0
        while j<donqty:
            mask.append(i)
            j+=1

    fn = 0
    start = Hbondlist(donliste[0][mask],Hliste[0],akzliste[0])
    Veranderung = [""]*len(donliste)
    vergleich1 = start

    while fn < len(donliste) :
        if fn%100 == 0:
            print(fn,end = '\r')
        vergleich2 = Hbondlist(donliste[fn][mask],Hliste[fn],akzliste[fn])
        k = 0
        while k < len(vergleich2):
            differenz = set(vergleich2[k]).symmetric_difference(set(vergleich1[k]))
            for d in differenz:
                Veranderung[k] += str(fn) +':' + str(d) + ';'
            k += 1
    
        vergleich1 = vergleich2
        fn += 1
    return Veranderung , start


def spatiallist(donliste,Hliste,akzliste,bins,radius,donqty=1):
    """
    Returns an array with shape (bins,frames) which contains the average number of Hydrogen bonds in the area around
    the z-axiswith equally distanced sections. The average number is then stored in the array with respect to the distance to the z-axis in the 
    middle of the box.
    Args:
        donliste: The subset of a Trajectory with the donor-atom of each donor-group.
        Hliste: The subset of a Trajectory with the hydrogen atom of each donor-group.
        akzliste: The subset of a Trajectory with all acceptor atoms.
        bins: Number of sections in which the the atoms are selected.
        radius: The maximal distance to the z-axis in the center of the box.
    """
    donorquantity=[]
    for i in range(len(donliste[0])):
        j = 0
        while j<donqty:
            donorquantity.append(i)
            j+=1
    frames = len(donliste)
    unit = np.array([0,0,1])
    fn = 0
    spatialarr = np.zeros((bins,frames))
    while fn < frames:
        if fn%100 == 0:
            print(fn,end = '\r')
        point = np.array([0.5*Hliste[fn].box[0][0],0.5*Hliste[fn].box[1][1],0])
        diff = Hliste[fn]-point
        distance=np.linalg.norm(diff-((diff)*unit)*unit,axis=1)
        areas = np.array(range(bins+1))/bins*radius
        masks =[]
        i = 0
        while i < bins:
            masks.append(((distance>=areas[i]) & (distance<areas[i+1])))
            i+=1
        liste = Hbondlist(donliste[fn][donorquantity],Hliste[fn],akzliste[fn])
        m = 0
        while m < bins:
            for k in liste[masks[m]]:
                spatialarr[m][fn]+=len(k)
            spatialarr[m][fn] = spatialarr[m][fn]/len(liste[masks[m]])
            m+=1
        fn+=1
    return spatialarr

def getHlifeplateus(frames,akzlist,sortedchangelist,mean=3):
    """
    Returns a boolean array with shape (akznr,frames). The array is the on and off time of hydrogen bonds a donor has 
    with each acceptor it has ever bonded with. If an on and off frequency is too high it can be averaged over a number
    frames. Also returns a list of frames in which the uneven frame denotes when a bond is formed and the even specify 
    when a bond is disbanded.(May contain smaller bugs)
    Args:
        frames: The number of frames in the Trajectory
        akzlist: The list of acceptors a specific donor has bonded with.
        sortedchangelist: The list of hydrogen bond changes for a specific donor sorted by the number of the acceptor
        mean: The minimal number of frames in which two hydrogen bonds are seperate.
    
    """
    sortedchange = sortedchangelist
    length = len(akzlist)
    frames = frames
    arr = np.zeros((length,frames),dtype=bool)
    k = 0
    Plateus = list()
    akzlen = len(akzlist)
    while k<akzlen:
        isBindung = False
        i = 0
        for t in sortedchange[k]:
            templis = [l.split(',') for l in t[0].split(':')]
            templisint = int(templis[0][0])
            while i< templisint:
                    arr[k][i] = isBindung
                    i+=1
            isBindung = not isBindung
        while i<frames:
            arr[k][i] = isBindung
            i+=1
        indices = np.where(arr[k] == True)
        n = 1
        indilen = len(indices[0])
        while n < indilen:
            if (indices[0][n] - indices[0][n-1])<mean:
                arr[k][indices[0][n-1]:indices[0][n]+1] = True
            n+=1
        l=1
        for sort in sortedchange:
            if len(sort)%2 != 0:
                del sort[-1]
        if arr[k][0] == True:
            Plateus.extend([0])
        if arr[k][frames-1] == True:
            Plateus.extend([frames])
        while l < frames:
            if (arr[k][l] != arr[k][l-1]):
                Plateus.extend([l])
            l+=1
        k+=1
    return Plateus,arr

def getaverageHbond(startlist,changelist,framenumber):
    """
    Return the an array with shape(donornr,framenumber) which contains how many hydrogen bonds each donor in each frame has.
    Args:
        startlist: The list of Hydrogenbonds of each donor-group for the first frame.
        changelist: The changes of hydrogen-bonds for each frame for each donor-group.
        framenumber: The number of frames in the trajectory.
    """
    lifetime = np.zeros(len(startlist),framenumber,dtype=np.int8)
    start = startlist
    k = 0
    while k<len(framenumber):
        sys.stdout.flush()
        if k%100 == 0:
            print(k,end='\r')
        Bindungen = set(start[k])
        temp =  [l.split(',') for l in changelist[k].split(';')]
        del temp[-1]
        i = 0
        for t in temp:
            temp2 = [l.split(',') for l in t[0].split(':')]
            while i< int(temp2[0][0]):
                lifetime[k][i] = len(Bindungen)
                i+=1
            if int(temp2[1][0]) in Bindungen:
                Bindungen.remove(int(temp2[1][0]))
            else:
                Bindungen.add(int(temp2[1][0]))
        while i<framenumber:
            lifetime[k][i] = len(Bindungen)
            i+=1
        k +=1
    return lifetime

def getlists(startlist,changelist,donnumber):
    """
    Returns the hydrogen bond changes of a specific donor sorted by the acceptor atom. Also returns a list with the number
    of all ever bonded acceptor atoms.
    Args:
        startlist: The list of Hydrogenbonds of each donor-group for the first frame.
        changelist: The changes of hydrogen-bonds for each frame for each donor-group.
        donnumber: The number of the donor which should be evaluated.
    
    """
    
    def sorting(key):
        kay = key[0]
        return int(kay[kay.find(':')+1:])
    
    temp =  [l.split(',') for l in changelist[donnumber].split(';')]
    if len(startlist[donnumber]) != 0:
        for i in startlist[donnumber]:
            temp.insert(0,['0:'+str(i)])
    del temp[-1]
    sor = sorted(temp, key=sorting)
    lis = list()
    liakz = list()
    i = 1
    k = 0
    while i< len(sor):
        if sor[i][0][sor[i][0].find(':')+1:] != sor[i-1][0][sor[i-1][0].find(':')+1:]:
            lis.append(sor[k:i])
            liakz.extend([sor[k][0][sor[k][0].find(':')+1:]])
            k = sor.index(sor[i])
        i+=1
    lis.append(sor[k:])
    liakz.extend([sor[k][0][sor[k][0].find(':')+1:]])
    return lis,liakz

def getHdistribution(Plateus):
    """
    Returns the lifetimes of hydrogen bonds of a donor.
    Args:
        Plateus: A list of on and off times of hydrogen bonds.
    """
    Times = np.array(Plateus[1::2]) - np.array(Plateus[::2])
    return Times

def helpfunction(i):
        sortedchangelists , akzlist = getlists(start,Veranderung,i)
        Plateus,arr = getHlifeplateus(frames,akzlist,sortedchangelists,mean= mean)
        return (getHdistribution(Plateus))
    
def parallelHdist(startlist,changelist,framenumber,average=3,cores=8):
    """
    Evaluates the hydrogen bond lifetime distributio for all frames and all donor-groups on multiple cpu-cores parallel.
    Args:
        startlist: The list of Hydrogenbonds of each donor-group for the first frame.
        changelist: The changes of hydrogen-bonds for each frame for each donor-group.
        framenumber: The number of frames in the trajectory.
        average: The minimal number of frames in which two hydrogen bonds are seperate.
        cores: The number of cpu-cores which should be used for evaluation
    """
    frames = framenumber
    start = startlist
    Veranderung = changelist
    mean = average
    final = collections.Counter()
    k = range(len(startlist))
    pool = multiprocessing.pool.Pool(cores)
    chunks = [k[x:x+120] for x in range(0, len(k), 120)]
    for c in chunks:
        print(c[0])
        results = pool.map(functools.partial(helpfunction,start=start,Veranderung=Veranderung,frames=frames,mean=mean),c)
        for r in results:
            final.update(r)
    pool.close()
    pool.join()
    return final
