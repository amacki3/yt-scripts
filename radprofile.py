import yt
import matplotlib.pyplot as plt
import filament
import numpy as np
import sys
import math
import profiling
import proj2
import matplotlib.patches as mpatches
#Load dataset

dataset = sys.argv[1]

ds = yt.load(dataset)

#Load filaments

filaments = filament.load_filaments(sys.argv[2])


print("Filaments read.")
##### For density change var to density and weight to cell_volume, remember to change graph labels!####
#####
#####

#####


#####



var="cylindrical_radial_velocity"


#for j in xrange(len(filaments)):
for j in xrange(2):
    print("looping")
#loop over all filaments
    profiles = []
    distlist = []
    used = []
    for i in xrange(1,filaments[j].shape[0]):
    #gather profile 


        profile = profiling.det_fil_profile(filaments[j][i-1],filaments[j][i],var,"cell_mass",2,ds)
    
        profiles.append(profile)
        used.append(profile.used)
    #determine the distance from start to this part of filament
        disttemp = []
        for s,e in zip(filaments[j][i-1],filaments[j][i]):
            if abs(s-e) > 0.5:
                disttemp.append( 1 + e - s)
            else:
                disttemp.append(s-e)
        disttemp = np.array(disttemp[:])  
        
        if len(distlist) > 0:
            distlist.append(distlist[-1] + np.sqrt(disttemp.dot(disttemp)))
        else:
            distlist.append(np.sqrt(disttemp.dot(disttemp)))
        
        del disttemp
                
    yarray = []
    number = 0
    distlist = ds.arr(distlist,'code_length')
    distlist = distlist.in_units('Mpc')

#loop over all profiles gathered and plot the variable against radius
    for profile in profiles:
        x = profile.x[profile.used]
        y = profile[(var)]
        yarray.append(y)
        y = y[profile.used]
        number = number + 1
        if number % 1 ==0:
            plt.plot(x,y,'k')
        
    plt.xlabel("Radius (Mpc)")
    plt.ylabel("Density (g/cm^3)")
    plt.xscale('log')
    #plt.yscale('log')
    plt.savefig("".join(["profiledens",str(j),".png"]))

    
    colors = ['b','r','g','k','m']
    handle = []
    for index in range(len(colors)):
        handle.append(mpatches.Patch(color=colors[index], label = "".join([str(profile.x[index*2].in_units('Mpc')), " Radius Profile"]))) 



    plt.clf()
    yavg = []
#now plot variation of variable along filament
    ind = 0
    for k in xrange(len(yarray[1])):
   # for k in xrange(1):
        y = []
        x = []


        for i in xrange(len(yarray)):
            if used[i][k]:
                x.append(distlist[i])
                y.append(yarray[i][k])
            if ind < 5:
                plt.plot(x,y,colors[k])
                ind = ind + 1
                yavg.append(y)
        yavgarray = np.array(yavg)
    del yavg,y
        
    plt.xlabel("Position on Filament (Mpc)")
    plt.yscale('log')
    plt.legend(handles=handle)
    if var == 'temperature':
        plt.ylabel("Temperature (K)")
        plt.fill_between(x,10E4,10E6,color='grey',alpha='0.2')
    else:
        plt.ylabel("Density (g/cm^3)")
    plt.savefig("".join(["tempdiff",str(j),".png"]))
    plt.clf()
    del x

    '''
    #now plot average of filament, use stddev as error
    yaverage = []
    yerror = []
    
    for index in range(len(profiles[0].x)):
        yaverage.append(np.mean(yavgarray[index]))
        yerror.append(np.std(yavgarray[index]))

    #plt.errorbar(profiles[0].x,yaverage,'k',yerr=yerror)
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel("Average Density of Filament (g/cm^3)")
    plt.xlabel("Distance (MPc)")
    plt.savefig("".join(["dmavg",str(j),".png"]))
    del yarray,distlist,y,x,profiles
    
    if j == 1:
        my_fils = filaments[1]
        plot = proj2.projplot(ds,filaments,'x')
        for i in xrange(len(filaments[1])):
            if i == 0:
                i = 1
            y1 = my_fils[i-1][1]
            y2 = my_fils[i][1]
            z1 = my_fils[i-1][2]
            z2 = my_fils[i][2]
        #Plot each individual segment
            if (abs(y2-y1) < 0.9) and (abs(z2-z1) < 0.9):
                plot.annotate_line((y1,z1),(y2,z2),coord_system='axis',plot_args={'color':'red'})
            del y1,y2,z1,z2
        plot.save("filshow.png")
    
        '''
    del distlist
