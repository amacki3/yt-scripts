from matplotlib import pyplot
import numpy as np
import sys
from scipy import stats
import yt
import os

def filament_length(filament):
    length = 0
    for i in range(filament.shape[0] - 1):
        p1 = filament[i][0:3]
        p2 = filament[i+1][0:3]
        length += np.sqrt(((p1 - p2)**2).sum())
    return length

def load_filaments(filename):
    start = filename.rfind("/",1) +1
    end = start + 6


    dsname = filename[start:end]
    fils = np.load(''.join(['/shome/mackie/data/',dsname,'/fils.npy']))
    return fils

def filter_filaments(filaments,keeplist):
    new_fils = [[] for i in range(len(filaments))]
    for i,fil in enumerate(filaments):
        for j,seg in enumerate(fil):
            if keeplist[i][j] == True:
                new_fils[i].append(seg)

    new_fils = np.array(new_fils)
    return new_fils


def load_keep_fils(filename):
    start = filename.rfind("/",1) +1
    end = start + 6

    dsname = filename[start,end]
    fils = np.load(''.join(['shome/mackie/data/',dsname,'/filkeep.npy']))
    return fils

def gen_filaments(filename):
    data = np.loadtxt(filename)

    p1 = data[:, 0:3]
    p2 = data[:, 3:6]
    f1 = data[:, 6]
    f2 = data[:, 7]

    d1 = np.rollaxis(np.concatenate([np.rollaxis(p1, 1), [f1]]), 1)
    d2 = np.rollaxis(np.concatenate([np.rollaxis(p2, 1), [f2]]), 1)
    del data, p1, p2, f1, f2

    all_filaments = []
    all_filaments.append([d1[0]])
    for i in range(1,d1.shape[0]):
        d1[i, 0:3] = np.mod(d1[i, 0:3], 1)
        d2[i-1, 0:3] = np.mod(d2[i-1, 0:3], 1)
        if (d1[i] == d2[i-1]).all():
            all_filaments[-1].append(d1[i])
        else:
            all_filaments[-1].append(d2[i-1])
            all_filaments.append([d1[i]])
    all_filaments[-1].append(d2[-1])

    print "Identified %d individual filaments." % len(all_filaments)
    all_filaments = [np.array(my_fil) for my_fil in all_filaments]

    keep = [True for my_fil in all_filaments]
    for i in range(len(all_filaments)):
        for j in range(i):
            check = all_filaments[i] == all_filaments[j]
            if isinstance(check, np.ndarray):
                if check.all():
                    keep[i] = False
            elif check:
                keep[i] = False            
            else:
                continue

    new_filaments = []
    for i, my_fil in enumerate(all_filaments):
        if keep[i]:
            new_filaments.append(my_fil)

    print "Identified %d unique filaments." % len(new_filaments)
    return new_filaments

if __name__ == "__main__":
    all_fil = load_filaments(sys.argv[1])
    segment_length = np.array([my_fil.shape[0] for my_fil in all_fil])
    all_length = np.array([filament_length(my_fil) for my_fil in all_fil])
	
    iseg = np.argmax(segment_length)


    with open("filout",'w') as file:
        for fils in all_fil:
           file.write("{}\n".format(fils))


    pf = load("/disk8/brs/core/double_core/full_box/sp_0_mackie/DD0046/DD0046")
	
