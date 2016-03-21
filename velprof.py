import os
import numpy as np
import matplotlib.pyplot as plt
import yt
import filament
import profiling
import filsep
import probmap
import yt.units as u
from yt import YTArray
from yt import YTQuantity
import parallel_analysis_interface as ytpar
yt.enable_parallelism()
ytpar.enable_parallelism()
def gen_dists(filament,ds):
    #Routine to determine the distances, in Mpc, between sections of a filament.
    
    distlist = np.empty(filament.shape[0])
    distlist[0] = 0.0
    
    for i in xrange(1,len(filament)-1):
        disttemp = (np.linalg.norm(filament[i] - filament[i -1]))
        if disttemp > 0.5:
            disttemp = 1.0 - disttemp 
        distlist[i] = disttemp + distlist[i-1]
    del disttemp
    #convert np array to yt dimensionful array, then convert to Mpc
    distlist = ds.arr(distlist,'code_length')
    distlist = distlist.in_units('Mpc')
    return distlist


def get_masses(dens_profiles,x,filaments,ds,accumulate=False):
    masses = YTArray(np.zeros((len(dens_profiles),10)),'g')
    stordict = {}
    fillen = np.arange(0,(len(filaments) -1 ))
    for sto,index in ytpar.parallel_objects(fillen,storage=stordict):
    #for index in range(len(filaments) - 1):
    #get cylinder object for all profiles
        normal,height,midpoint = profiling.get_cyl_params(filaments[index],filaments[index+1])
        
        height = height * ds.length_unit


        #We check to see if the mass should be accumulated, or if it ought to be reported only in shells
        #Then we find the correct volume relating to that.

        shell = np.roll(x,1)
        shell[0] = 0
        vol = np.pi * height * ( (x**2) - (shell**2))
        

        mass = ( (vol  * dens_profiles[index] ).in_units('g'))

        if accumulate == True:
            for i in range(1,mass.size):
                mass[i] = mass[i] + mass[i-1]

        
        sto.result = mass
        
        sto.result_id = index
        
        del vol,mass

    #Recombine results from parallel processing
    for key,values in sorted(stordict.items()):
        masses[key] = values
        

    return masses



def plot_vel(filament,ds,dataset,fil=-1,maskarray=False):
    #Routine to plot velocity of particles along a filament
    #Set gravitational constant
    G = YTQuantity(6.67408E-11,'m**3/(kg * s**2)')
    #Gather velocity and density values from disk, done in parallel to speed computation
    #We first gather a list of the profiles on the disk, then reshape this into a list of [density,other] profiles
    #This is then iterated over in parallel to load the correct data
    filelist = sorted(os.listdir(''.join(['/shome/mackie/data/',dataset,'/profiles'])))
    profnumbers = len(filelist)/2
    files = [ [filelist[i],filelist[i+profnumbers]] for i in range(profnumbers)] 
    del filelist,profnumbers
    



    storage = {}
    for stor,file_in_dir in ytpar.parallel_objects(files,storage=storage):
        
        #Determines correct index to give to segment profiles
        filnum = int(file_in_dir[0][7:10])
        segnum = int(file_in_dir[0][13:16])
        #Calc total density for each segment 
        densprof = yt.load(''.join(['/shome/mackie/data/',dataset,'/profiles/',file_in_dir[0]]))
        dm = densprof.data['dark_matter_density'].in_units('g/cm**3')
        dens = densprof.data['density'].in_units('g/cm**3')
        totaldens = dm + dens
        del densprof,dm,dens
        #Get velocity profiles
        velprof = yt.load(''.join(['/shome/mackie/data/',dataset,'/profiles/',file_in_dir[1]]))
        vel = velprof.data['cylindrical_radial_velocity'].in_units('km/s')

        stor.result = (vel,totaldens)
        stor.result_id = (filnum,segnum)



        
    #Restruct dict into np array of correct structure.
    vel_profs = [ [] for i in range(len(filament))]
    densprofs = [ [] for i in range(len(filament))]
    x = yt.load(''.join(['/shome/mackie/data/',dataset,'/profiles/',file_in_dir[0]])).data['x']
    xarr = [[] for i in range(len(filament))]
    for key,values in sorted(storage.items()):
        filnum,segnum = key
        vel,dens = values
        xarr[filnum].append(x.in_units('Mpc'))
        
        vel_profs[filnum].append(vel.in_units('km/s'))
        densprofs[filnum].append(dens)
    for i in range(len(xarr)):
        xarr[i] = YTArray(np.array(xarr[i]),'Mpc')
        vel_profs[i] = YTArray(np.array(vel_profs[i]),'km/s')
    vel_profs = YTArray(np.array(vel_profs),'km/s')
    xarr = YTArray(np.array(xarr),'Mpc')
    del storage
    #Turn into np arrays for QoL
    densprofs = np.array(densprofs)
    #Gather x bins from disk
    
    #Determine masses and thus escape velocities
    
      
    mass = [get_masses(densprofs[i],x,filament[i],ds,accumulate=True) for i in range(len(filament))]

    mass = np.array(mass)
    mass = YTArray(mass,'g')
    print mass[1][1]
    
    del densprofs
    
    print xarr[1][1]

    vel_ratio = ( ( (2*G* mass) / xarr) ** (1.0/2.0))
    vel_ratio = vel_ratio.in_units('km/s')

    if yt.is_root():                        
        print mass[1][1]
        print xarr[1][1]
        print vel_ratio[1][1]
        
#vel_ratio is **approx** escape vel, used to ratio later
    #Generate ratio of velocity to escape velocity
    vel_profs = (vel_profs.in_units('km/s')/vel_ratio.in_units('km/s'))
    del vel_ratio
    



    if fil > -1:
        
        length_plot =  plot_vel_fil(vel_profs[fil].v,gen_dists(filament[fil],ds),x)
    else:
        length_plot = None

    if maskarray:
        print"Masking Data"
        vel_profs = np.ma.masked(vel_profs, mask = ~maskarray,fill_value=np.nan)
    #Flatten vel profs, ought to bemore elegant solution
    vel_prof_flatten = []    
    for fil in vel_profs:
        for seg in fil:
            vel_prof_flatten.append(seg)
    vel_profs = np.array(vel_prof_flatten)
    del vel_prof_flatten

    plot = probmap.prob_heat_map(vel_profs,'radial_velocity',x=x)
    
    return plot,length_plot



def plot_vel_fil(vel_prof,dists,y):
    #Generates plot of velocities along a filament
    #Generate plot
    if yt.is_root():
        plot = plt.figure()
        #Set up x and y values, and correctly fill and ensure shape is correct
        x = np.empty((dists.shape[0],y.shape[0]))
        yy = np.empty((vel_prof.shape[0],y.shape[0]))
        x = x.T
        for i in xrange(vel_prof.shape[0]):
            yy[i] = y
        for i in xrange(y.shape[0]):
            x[i] = dists
        x = x.T
        #Ensure our colorbar is in correct form
        colormax = np.nanmax(vel_prof)
        colormin = abs(np.nanmin(vel_prof))
        colormax = max([colormax,colormin])
        colormin = 0 - colormax

        print( x.shape, yy.shape, vel_prof.shape)
                
        #Create a 'probability' map plot
        vel_map = plt.pcolormesh(x,yy,vel_prof,cmap='seismic', vmin =colormin, vmax=colormax)
        #Label graph
        plt.xlabel("Distance along filament (Mpc)")
        plt.ylabel("Radius from Filament (Mpc)")
        cbar = plt.colorbar(vel_map)
        cbar.ax.set_ylabel("Velocity (km/s)")
  
    
    
        return plot
    
if __name__ == "__main__":
    import sys
    
    filaments = filament.load_filaments(sys.argv[1])
    ds = yt.load(sys.argv[1])
    dataset = sys.argv[1][sys.argv[1].rfind("/",1):]

    fil = 1
    
    dataset=dataset[1:]
    print dataset
    
    plot,lplot = plot_vel(filaments,ds,dataset)
    plot.savefig("".join(["velplot",dataset,".png"]))
    if lplot is not None:
        lplot.savefig("".join([dataset,'velalong',fil,".png"]))
