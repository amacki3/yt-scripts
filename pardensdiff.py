import sys
import os
import numpy as np
import yt
import filament
import probmap
import profiling
import filsep
yt.enable_parallelism()

def gen_profs(ds,fils,dsname,keep=True,totalratio=None,mask=False):
    x = np.zeros(10)
    

    
    storage = {}
    #Check to see if totalratio has been set, else determine it.
    if totalratio == None:
    #Determine simulation-wide ratio of baryon to Total matter, done in parallel to speed up
        for stor,i in yt.parallel_objects(range(1),storage = storage):
                        
            denstotal = ds.all_data().quantities.weighted_average_quantity("density","cell_volume")
            dmtotal = ds.all_data().quantities.weighted_average_quantity("dark_matter_density","cell_volume")
            
            totalratio = denstotal / (dmtotal + denstotal)
            del dmtotal

            stor.result = totalratio,denstotal.in_units('g/cm**3')
            stor.result_id = 1

        totalratio,denstotal = storage[1]
    print totalratio # print to stdout, allows us to check ratio makes sense/debug
    del storage
    #Print to a file, allows checking of progress when in job queue
    if yt.is_root():
        with open("Testinfo.dat",'w') as f:
            f.write('Ratio determined')
    print("All Ratio determined")
    x = []

    storage = {}
    #Gather list of all profiles on disk
    filelist = sorted(os.listdir(''.join(['/shome/mackie/data',dsname,'/profiles'])))
    #Include only density profiles
    filelist = filelist[:len(filelist)/2]
    #Load a numpy mask array if required
    if keep == True:
        keeplist = np.load(''.join(['/shome/mackie/data',dsname,'/filkeep.npy']))
    
    
    for stor,file_in_dir in yt.parallel_objects(filelist,storage=storage):
        
        filnum = int(file_in_dir[7:10])
        segnum = int(file_in_dir[13:16])

        prof = yt.load(''.join(['/shome/mackie/data',dsname,'/profiles/',file_in_dir]))
        dm = prof.data['dark_matter_density']
        dens = prof.data['density'].in_units('g/cm**3')
        
        result = (dens /( dm + dens)) - totalratio
        dens = dens/denstotal
        result = result.v
        if keep == True:
            if keeplist[filnum][segnum] == True:    
                stor.result = result,dens
                stor.result_id = file_in_dir
        else:
            stor.result = result,dens
            stor.result_id = file_in_dir
        del prof,dm,dens,filnum,segnum
    
    results = []
    denresult = []
    for keys,values in storage.items():
        results.append(values[0])
        denresult.append(values[1])
    results = np.array(results)
    denresult = np.array(denresult)

    if mask:
        results = np.ma.masked_array(results,mask=~mask,fill_value=np.nan)
#Get x bins
    x = yt.load(''.join(['/shome/mackie/data',dsname,'/profiles/densfil000seg000.h5'])).data['x']

    if yt.is_root():
        with open("Testinfo.dat",'a') as f:
            f.write("returning results to plot")

        print results
    
    return results,denresult,x
    




'''



    storage = {}
    






    for results,(i,filaments) in yt.parallel_objects(enumerate(fils),storage=storage):
        
        profiles = []

        for seg_id  in xrange(1,filaments.shape[0]):
            
            profile = profiling.det_fil_profile(filaments[seg_id - 1],filaments[seg_id],[("dark_matter_density"),("density")],"cell_volume",2,ds)
            
            dmprof = np.ma.masked_values(profile[("dark_matter_density")].v,value=0.0,atol=1e-40)
            baryon = np.ma.masked_values(profile[("density")].v,value=0.0,atol=1e-40)

            result = (baryon / (dmprof + baryon) ) - totalratio
            profiles.append(result)
            
            x = profile.x
        
        results.result = profiles
        results.result_id = i

        print 'filament number:'
        print i

        print ("Results gen. Now plotting")

       

    results = [ result for (result,key) in sorted(results.items())]
   
    return results,x
            
'''            


if __name__ == "__main__" :
    
    #Load Dataset

    dataset = sys.argv[1]
    ds = yt.load(dataset)
    filaments = filament.load_filaments(sys.argv[1])
    dsname = dataset[dataset.rfind('/',1):]

    profiles,dens,x = gen_profs(ds,filaments,dsname,keep=False)
    #need to flatten profiles

    if yt.is_root():
        with open('info.dat','w') as f:
            f.write("Profiles generated. Now gening prob map")
        
        prob = probmap.prob_heat_map(profiles,'densdiff',x)
        prob.savefig(''.join([str(ds.current_redshift),'zdendiff.png']))
        prob = probmap.prob_heat_map(dens,'densrat',x)
        prob.savefig(''.join([str(ds.current_redshift),'zdens.png']))


        with open("densitydiffs.dat",'w') as f:
            for i in range(len(profiles)):
                for j in range(len(profiles[i])):
                    if abs(profiles[i][j]) > 0.05:
                        f.write(''.join([str(i),"/t",str(j),"/n"]))

