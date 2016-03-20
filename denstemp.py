import sys
import parallel_analysis_interface as ytpar
import yt
import matplotlib.pyplot as plt
import filament
import numpy as np
import velprof
import os
from yt import YTArray
from mpi4py import MPI


ytpar.enable_parallelism()
yt.enable_parallelism()

def plot_prob2D(profiles1,profiles2,temp=False):
    
    #Function which takes two lists of profiles and plots a probability map for the two variables

    min1 = np.nanmin(np.ma.masked_equal(profiles1,0))
    min2 = np.nanmin(np.ma.masked_equal(profiles2,0))
    max1 = np.nanmax(profiles1)
    max2 = np.nanmax(profiles2)
    print profiles1.size
    print min1,max1
    print min2,max2

    #Set up bins for plot
    bins = 50
    x = np.logspace(np.log10(min1),np.log10(max1),bins)
    y = np.logspace(np.log10(min2),np.log10(max2),bins)
    
    #Set up array of values to populate probability map
    prob = np.zeros([x.shape[0],y.shape[0]])

    
    print(x,y)
    #Determine where values should reside in probability array
    #Requires looping through all elements in both profile matrices and then allocating an increment to correct area of prob array.

    for i in range(profiles1.shape[0]):
        for j in range(profiles1[i].shape[0]):
            x_index = -1
            y_index = -1
            if not(profiles1[i][j] == np.nan) and not(profiles2[i][j] == np.nan):
                for index in range(x.size):
                    if x_index < 0 and profiles1[i][j] < x[index]:
                        x_index = index
                    if y_index < 0 and profiles2[i][j] < y[index]:
                        y_index = index
                        break

            prob[x_index][y_index] += 1

    print(prob)
    prob = prob / np.sum(prob)

    plot = plt.figure()
    CS = plt.pcolormesh(x,y,prob.T,cmap = 'YlOrRd')



    plt.xlabel("Density (g/cm^3)")
    plt.ylabel("Temperature (K)")

    plt.xscale('log')
    plt.yscale('log')

    plt.fill_between(x,10E4,10E6,color='grey',alpha='0.2')

    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel("Probability")
    return plot
            
if __name__ == "__main__":

    #Load filaments and dataset
    filaments = filament.load_filaments(sys.argv[2])
    ds = yt.load(sys.argv[1])
    dsname = sys.argv[1][sys.argv[1].rfind("/",1):]

    #Get profiles for density and temperature
    #mask = np.load('/shome/mackie/data/DD0252/tempkeep.npy')
   


    #Get profiles for density and temperature
    filelist = sorted(os.listdir(''.join(['/shome/mackie/data/',dsname,'/profiles'])))
    profnumbers = len(filelist)/2
    files = [ [filelist[i],filelist[i+profnumbers]] for i in range(profnumbers)]
    #Filter out 'bad' filament segments
    keeplist = list(np.load(''.join(['/shome/mackie/data',dsname,'/filkeep.npy'])))   
    
    filaments = filament.filter_filaments(filaments,keeplist)


    #need to gather density and temp profiles
    x = yt.load(''.join(['/shome/mackie/data',dsname,'/profiles/',filelist[0] ] )).data['x']
    del filelist,profnumbers

    storagedict = {}
    
    
    for stor,file_in_dir in ytpar.parallel_objects(files,storage=storagedict):
        
        filnum = int(str(file_in_dir[0][7:10]))
        segnum = int(str(file_in_dir[0][13:16]))
        velfil = file_in_dir[0]
        

        
        prof = yt.load(''.join(['/shome/mackie/data',dsname,'/profiles/',file_in_dir[0]]))
        
        dens = prof.data['density'].in_units('g/cm**3')
        
        temppro = yt.load(''.join(['/shome/mackie/data',dsname,'/profiles/',file_in_dir[1]]))
        temp = temppro.data['temperature'].in_units('K')
        
        
        if keeplist[filnum][segnum]:
            stor.result = dens,temp
            stor.result_id = velfil
        
        del filnum,segnum,velfil,prof,dens,temppro,temp

    print 'COMPLETED LOOP'
    temp = [[] for i in range(len(filaments))]
    temporarydens = [[] for i in range(len(filaments))]

    for keys,values in sorted(storagedict.items()):
        if not values ==  None:
            
            dens,temperature = values
            filnum = int(str(keys[7:10]))
        
            temp[filnum].append(YTArray(temperature,'K'))
            temporarydens[filnum].append(YTArray(dens,'g/cm**3'))
    if yt.is_root():
            print dens

    del storagedict
    totaldensity = 0
    outtemp = 0 

    #mass = velprof.get_masses(temporarydens[0],x,filaments,ds)
    
    mass = [velprof.get_masses(temporarydens[i],x,fil,ds) for i,fil in enumerate(filaments)]
    print 'mass generated'
    dens = []
    tempflat = []
    if np.any(mask):
        mass = np.ma.masked_where(mask==False,mass)
        np.ma.set_fill_value(mass,0.0)

    tempkeep = [[] for i in range(len(filaments))]
    for i,fil in enumerate(temp):
        for j,seg in enumerate(fil):
            tempkeep[i].append([ True for x in range(10) ])
            dens.append(temporarydens[i][j])
            tempflat.append(temp[i][j])

            for k,t in enumerate(seg):
                if 1E5< t < 1E7:
                    if np.any(mask):
                        if mask[i][j][k] == False:
                            dens[-1][k] = np.nan
                            tempflat[-1][k] = np.nan
                    outtemp = mass[i][j][k] + outtemp
                    tempkeep[i][j][k] = False
                totaldensity = totaldensity + mass[i][j][k]
    tempkeep = np.array(tempkeep)
    del temporarydens
            
    if yt.is_root():
        
        with open('/shome/mackie/TEMPERATUREDETAILS.txt','a') as f:
            f.write('z: %f \t ratio: %f \n' %(ds.current_redshift,outtemp/totaldensity))
        print(outtemp/totaldensity)
        print('above ratio of invis matter to all matter, below tot. of invis and then total')
        print(outtemp)
        print(totaldensity)

        plot = plot_prob2D(np.array(dens),np.ma.masked_equal(tempflat,0.0))
        plot.savefig(''.join([str(ds.current_redshift),'tempdensmap.png']))
        
        if not mask:
            tempkeep = np.array(tempkeep)
            np.save(''.join(['/shome/mackie/data/',dsname,'/tempkeep.npy']),tempkeep)


    if yt.is_root():
        print("Root checking does actually work")
