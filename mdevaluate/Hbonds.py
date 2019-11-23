import numpy as np
import scipy
from mdevaluate.pbc import pbc_diff
import sys
import collections
from functools import lru_cache

def recursive_len(item):
        if type(item)==str:
            return len(item)
        if ((type(item) == list) or (type(item) == tuple) or (type(item)==dict)):
            return sum(recursive_len(subitem) for subitem in item)
        else:
            return 1
class MyList(list):
    def __init__(self, *args):
        list.__init__(self, *args)
    def __repr__(self):
        return "Number of charakters = "+str(recursive_len(list(self)))

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
        donqty: The number of hydrogen each donor atom has.
        Returns during the calculation the current frame in steps of 100.
    """
    mask=[]
    for i in range(len(donliste[0])):
        j = 0
        while j<donqty:
            mask.append(i)
            j+=1
    fn = 0
    start = Hbondlist(donliste[0][mask],Hliste[0],akzliste[0])
    Veranderung = [""]*len(donliste[0][mask])
    vergleich1 = start
    while fn < 1000 :
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
    return MyList(Veranderung) , start


def spatiallist(donliste,Hliste,akzliste,bins,origin,radius,donqty=1,maxheight=np.Inf,):
    """
    Returns an array with shape (bins,frames) which contains the average number of Hydrogen bonds in the area around
    the z-axis with equally distanced sections. The average number is then stored in the array with respect to the distance to the z-axis in the 
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
        point = origin
        diff = Hliste[fn]-point
        distance=np.linalg.norm(diff-((diff)*unit)*unit,axis=1)
        areas = np.array(range(bins+1))/bins*radius
        masks =[]
        i = 0
        while i < bins:
            masks.append(((distance>=areas[i]) & (distance<areas[i+1]) & (Hliste[fn][:,2]<maxheight)))
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

def getaverageHbond(startlist,changelist,numberofframes):
    """
    Return the an array with shape(donornr,framenumber) which contains how many hydrogen bonds each donor in each frame has.
    Args:
        startlist: The list of Hydrogenbonds of each donor-group for the first frame.
        changelist: The changes of hydrogen-bonds for each frame for each donor-group.
        numberofframes: The number of frames in the trajectory.
    """
    lifetime = np.zeros((len(startlist),numberofframes),dtype=np.int8)
    start = startlist
    k = 0
    while k<len(startlist):
        sys.stdout.flush()
        if k%100 == 0:
            print(k,end='\r')
        Bindungen = set(start[k])
        temp =  [l.split(',') for l in changelist[k].split(';')]
        del temp[-1]
        i = 0
        for t in temp:
            temp2 = [l.split(',') for l in t[0].split(':')]
            temp2 = [l.split(',') for l in t[0].split(':')]
            while i< int(temp2[0][0]):
                lifetime[k][i] = len(Bindungen)
                i+=1
            if int(temp2[1][0]) in Bindungen:
                Bindungen.remove(int(temp2[1][0]))
            else:
                Bindungen.add(int(temp2[1][0]))
        while i<numberofframes:
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

def HlifePlateus(startlist,changelist,frames,average=3):
    """
    Returns the lifetime of all Hbonds of the system in a Counter object.
    Args:
        startlist: The list of Hydrogenbonds of each donor-group for the first frame.
        changelist: The changes of hydrogen-bonds for each frame for each donor-group.
        frames: The number of frames in your Trajectory.
        average: The minimum number of frames two Hbonds must be seperated to count as seperate Hbonds.

    """
    start=startlist
    Veranderung=changelist
    frames=frames
    mean = average
    final = collections.Counter()
    j=0
    while j<frames:  
        if j%100==0:
            print(j,end='\r')
        sortlist, akzlist = getlists(start,Veranderung,j)
        sortedchange = sortlist
        length = len(akzlist)
        arr = np.zeros((length,frames),dtype=bool)
        k = 0
        Plateus = list()
        akzlen = len(akzlist)
        while k<akzlen:
            isBindung = False
            i = 0
            for t in sortedchange[k]:
                templisint = int(t[0][:t[0].index(':')])
                arr[k][i:templisint] = isBindung
                i = templisint
                isBindung = not isBindung
            arr[k][i:frames] = isBindung
            indices = np.where(arr[k] == True)
            n = 1
            indilen = len(indices[0])
            while n < indilen:
                if (indices[0][n] - indices[0][n-1])<mean:
                    arr[k][indices[0][n-1]:indices[0][n]+1] = True
                n+=1
            l=1
            Plateus.append(np.where(np.diff(arr[k]))[0]+1)
            if arr[k][0]:
                Plateus[k] = np.insert(Plateus[k],0,0,axis=0)
            if arr[k][frames-1] == True:
                Plateus[k] = np.append(Plateus[k],[frames])
            k+=1
        for m in Plateus:
            final.update(m[1::2]-m[::2])
        j+=1
    return final

class Hbond():
    
    def __init__(self, donlist, Hlist, akzlist,donqty=1):
        self.akz=akzlist
        self.don=donlist
        self.H=Hlist
        self.donqty=donqty
        self.frames=len(Hlist)
        self.change,self.start=Veranderungsliste(Hliste=Hlist,akzliste=akzlist,donliste=donlist,donqty=donqty)
    @lru_cache(maxsize=10)
    def averageHbond(self):
        """
        Returns the average number of hydrogen bonds per donor group over the whole trajectory for each frame.
        """
        return np.mean(getaverageHbond(changelist=self.change,startlist=self.start,numberofframes=self.frames),axis=1)
    @lru_cache(maxsize=10)
    def averageHbondmol(self):
        """
        Returns the average number of hydrogen bonds for each donor group over the whole trajectory for each frame.
        """
        return getaverageHbond(changelist=self.change,startlist=self.start,numberofframes=self.frames)
    @lru_cache(maxsize=10)
    def lifetimedist(self,average=3):
        """
        Returns a Counter object with the number of lifetimes of the hydrogen bonds in the trajectory.
        Args:
            average: The minimum frame difference between a hydrogen bonding to count it as seperate.
        """
        return HlifePlateus(startlist=self.start,changelist=self.change,average=average,frames=self.frames)
    @lru_cache(maxsize=10)
    def Hbondlistmol(self,donnumber):
        """
        Returns two lists of a specific donor group. The first list returns all changes of the donor ordered by acceptor numbers,
        the second list are the numbers of all acceptors which have been interacted with.
        Args:
            donnumber: The number of the donor group which lists are calculated.
        """
        return getlists(changelist=self.change,startlist=self.start,donnumber=donnumber)
    
    
