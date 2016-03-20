import numpy as np
import matplotlib.pyplot as plt
import math
import yt


def prob_heat_map(profiles,var,x=np.empty(10)):
    #routine to generate a probability heat map for a variable(e.g. density/temperature) in a filament. Will also plot mean and median to this map.
    
    '''
    if 'densdiff' not in var:
    #set up x values, and create empty arrays
        x = profiles[0].x
        data = np.empty([len(profiles),x.shape[0]])
        mask = np.ones([len(profiles),x.shape[0]],dtype=bool)
    
    
    #populate data, and generate a mask for 0 values
        for index in xrange(len(profiles)):
            data[index,] = profiles[index][(var)].v
    
    else:
    '''
    data = np.array(profiles)
    if not("velocity" in var or "densdiff" in var or "densrat" in var):
        data = np.ma.masked_values(data,0.0,atol = 1e-40)
     
    
    data=np.ma.masked_invalid(data)

    
#determine min, max and the binsize for variable
    if yt.is_root(): print(np.ravel(data))
    mind = np.nanmin(data.reshape(-1))
    maxd = np.nanmax(data.reshape(-1))
    print 'max'
    print maxd
    #if mind < 0:
    #   mind = 0 - np.amax([abs(mind),maxd])
    #    maxd = 0 - mind
    bins = 50
    length = maxd - mind
#generate bins for the variable
    if var == "cylindrical_radial_velocity" or var == "densdiff" or var == "densrat":
        y = np.linspace(mind,maxd,bins)      
    else:    
        print mind
        if mind == 0:
            mind = data.min()
        print mind
        y = np.logspace(np.log10(mind),np.log10(maxd),bins)
   #set up prob array
    prob = np.zeros([10,y.size])
#loop over all values of variable and add to relevant bin in prob array

    for i in range(len(data)):
        for j in range(10):
            z = data[i][j]
            if math.isnan(z):
                break
            yindex = 0
    #calculate correct index to increment
            for index in range(y.size):
                if z <= y[index]:
                    yindex = index
                    prob[j,yindex] += 1
                    break
                
    
    mean = np.empty(x.size)
    median = np.empty(x.size)
    #normalize probability array, calculate mean and median
    for index in range(x.size):
        yvals = np.array(prob[index,:])
        summed = np.sum(yvals)
        if summed > 0:
            yvals = yvals / summed
        
        
            prob[index,:] = yvals
        del yvals,summed

    
    for index in range(x.size):
        maskeddata = data[:,index]
        mean[index] = np.nanmean(maskeddata)
        median[index] = np.nanmedian(maskeddata)
    del maskeddata
        
        
        
    
    #now plot
    plot = plt.figure()
    plt.plot(x,mean,'g',label='Mean')
    plt.plot(x,median,'r',label='Median')
   
    CS = plt.pcolormesh(x,y,prob.T,cmap='seismic')
    
    plt.xlabel("Radius (MPC)")
    plt.xscale('log')
    #plt.yscale('log')

    if var.lower() == 'density':
        plt.ylabel("Density (g/cm^3)")
        plt.yscale('log')
        plt.xlim(-0.1,2.5)
        

    elif var.lower() == 'temperature':
        plt.ylabel("Temperature (K)")
        plt.fill_between(x, 10E4 , 10E6 ,color='grey',alpha='0.2')
        plt.yscale('log')
        plt.xlim(-0.1,2.5)
            
    elif var.lower() == 'cylindrical_radial_velocity':
        plt.ylabel("Velocity/Escape Velocity")
        plt.xscale('linear')
        plt.yscale('linear')
        
    elif var.lower() == "densrat":
        plt.ylabel("Local Density/Average Density")
        plt.yscale('linear')
        plt.xlim(-0.1,2.5)

    else:
        plt.ylabel("Baryon Density/ Total Density - Avg. Baryon Fraction")
        plt.yscale('linear')
        plt.xlim(-0.1,2.5)

    
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel("Probability")
    plt.legend()
    return plot



if __name__== "__main__":
    import filament
    import sys
    import profiling

    filaments = filament.load_filaments(sys.argv[2])
    ds = yt.load(sys.argv[1])
    print("Filaments loaded.")
    
    vars = ['density','temperature','cylindrical_radial_velocity']
    
    #for var in vars:
    var = 'cylindrical_radial_velocity'
    if var == 'cylindrical_radial_velocity':
       # for fil in range(len(filaments):
        for fil in range(2):   
            profiles = []
    #for var in vars:
            for i in xrange(1,filaments[fil].shape[0]):
                profile = profiling.det_fil_profile(filaments[fil][i-1],filaments[fil][i],var,"cell_mass",2,ds)
                profiles.append(profile)
            print("profiling complete")
            plot = prob_heat_map(profiles,var)
            plot.savefig("".join(["prob",var[:5],str(fil),"mapt.png"]))

