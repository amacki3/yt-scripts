import sys
import os
import numpy as np
import yt
import filament
import profiling 
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from yt import YTArray
yt.enable_parallelism()



def generate_halos(ds,dsname):
    #function to find the halos for the given dataset, using the HOP finder
    #Calls a fn to gather positions and radii of halos over a mass threshold
    #Autosaves the HC
    hc = HaloCatalog(data_ds=ds,finder_method='hop',finder_kwargs={'threshold':300,'padding':0.04},output_dir="".join(['~/',dsname]))
    
    hc.add_filter('quantity_value','particle_mass','>',1E12,'Msun')
    hc.create()
    hc = hc.halos_ds.all_data()
    

    radii, positions = get_halo_pos(ds,hc)
    return hc, radii, positions

def load_halos(ds,location):
    #Merely loads an already saved halocatalog
    #Calls a fn to gather positions and radii of halos over a mass threshold
    halo_ds = yt.load(location)
    hc = HaloCatalog(data_ds=ds,halos_ds=halo_ds)
    hc.add_filter('quantity_value','particle_mass','>',1E12,'Msun')
    halo_ds = hc.halos_ds.all_data()

    radii, positions = get_halo_pos(ds,halo_ds)
    return hc, radii, positions
    



def get_halo_pos(ds,hc):
    
    halo_ds = hc

    radii = halo_ds['virial_radius']
    radii = radii / ds.length_unit
    radii = radii.v
    
    #Generates useful format of the positon and radius of each halo
    temppos = halo_ds[('particle_position_x')] / ds.length_unit
    temppos = temppos.v
    positions = np.empty((3,temppos.size))
    positions[0] = temppos
    temppos = halo_ds[('particle_position_y')] / ds.length_unit
    positions[1] = temppos.v
    temppos = halo_ds[('particle_position_z')] / ds.length_unit
    positions[2] = temppos.v
    del temppos
    
    return radii, positions.T


def sep_fils(filament,radii,positions,ds):
    #Routine which checks to see if segments of a filament reside within a halo



    #Init vars
    storage = {}
    fils = np.array(filaments)
    keep_arr = []

    #Loop across all filaments, utilizing parellization, split into pool of 4 jobs, this is arbritrary
    for sto,index in yt.parallel_objects(range(len(fils)),storage = storage,njobs=4):
        

        #init some new vars, keeplist will be a 'mask' array for filaments
        keeplist = np.empty(fils[index].shape[0],dtype=bool)
        segstor = {}
        #loop over all segments of a filament, again utlize parellization.
        for stor,(i,section) in yt.parallel_objects(enumerate(fils[index]),storage=segstor):
            
            keep = True
            
            #Check if segment is within halo or not
            for pos,rad in zip(positions,radii):

                if np.linalg.norm(section[0:3] - pos) < rad:
                    keep = False
                    break
            
            stor.result = keep
            stor.result_id = i

        for i in range(fils[index].shape[0]):
            keeplist[i] = segstor[i]


      

        sto.result = keeplist
        sto.result_id = index

    #Aggragate arrays from all parallel jobs
    keep_arr = [ lists for (key,lists) in sorted(storage.items())]
                          
    return keep_arr


def make_save_profiles(fils,dsname,ds,check=False):
    
#Routine which generates temperature, density and velocity profiles for all segments of a filament
    
    #inits useful vars
    
    var = [('density'),('dark_matter_density'),('temperature'),('cylindrical_radial_velocity')]
    weights = ("cell_volume","cell_mass")
    profiles = [[] for i in range(len(fils))]
    filsize = len(fils)


    #check to see if some profiles already exist on disk, in whcih case we dont need to make the profiles
    
    if check == True:
        checkfils = []
        directory = ''.join(['~/data',dsname.upper(),'/profiles/'])
        filelist = os.listdir(directory)
        for files in filelist:
            checkfils.append((int(files[7:10]),int(files[13:16])))
    else:
        checkfils = []

    #Loops over all filaments
    for i in yt.parallel_objects(range(filsize),dynamic=check):
        profs = []
        filaments = fils[i]

        for section in (range(len(filaments) - 1)):
            #Generate useful profiles for all segments of filament
            if (i,section) in checkfils:
                var_prof = profiling.det_fil_profile(filaments[section],filaments[section + 1],var[0:2],weights[0],2,ds)
                var_prof2 = profiling.det_fil_profile(filaments[section],filaments[section + 1],var[2:],weights[1],2,ds)
            

                #Save these to disk
                save_prof(var_prof,dsname,'dens',(i,section))
                save_prof(var_prof2,dsname,'kine',(i,section))


            #profs.append([var_prof,var_prof2])
            

'''
ERRORS IN PICKLING MAKE THIS LOCAL STORAGE UNABLE TO WORK, WILL LOOK INTO WORKAROUNDS
        stor.result = var_prof
        stor.result_id = i
    #Aggragate results
    profiles = [values for (key,values) in sorted(stordict.items())]
    
          
    return profiles
'''


def save_prof(profile,dsname,var,seg_id):
    #Saves a single profile
    #Requires a directory name(dsname), the var name to save under, a tuple of two ids, a filament and segment identifier. 
    import os        

    #Determine directory to print to
    directory = ''.join(['~/data',dsname.upper(),'/profiles/'])
    #Check dir exists,if not create dir
    try:
        os.makedirs(directory)
    except OSError:                     
        pass

    directory = ''.join([directory,var,'fil',str(seg_id[0]).zfill(3)])
    profile.save_as_dataset(''.join([directory,'seg',str(seg_id[1]).zfill(3)]))



def save_all_profs(profiles,dsname,ds):
    #Save ALL profiles to disk, where profiles are supplied in a list of 2 variable types ( those weighted by mass and those by volume), containing lists of filaments containing the list of segments.
    
    #Init useful vars, important varshort is 4 chars long for loading from disk later
    varshort = ['dens','kine']
        
    #Split job into 2, one for the different weight fields
    for jobnumber in yt.parallel_objects(range(2),njobs=2):
        #Parallelize the workload efficiently
        for i,fil in yt.parallel_objects(enumerate(profiles)):
         
            for j,seg in enumerate(fil):
               #Determine which variable we are dealing with
                if jobnumber % 2 == 0:
                    segment = seg[0]
                else:
                    segment = seg[1]

                #Save segment of filament profile to disk
                save_prof(segment,dsname,varshort[jobnumber],(i,j))
                
def load_all_profs(directory,filament_num):
    #Loads profiles from disk.
    #Profiles are stored as a 3D list of segments listed within filaments listen in var type
    #This means there are three keys needed to identify one profile, its var type, its filament number and its segment number
    #I.E. They are accessed via profiles[ variable number ] [ filament number] [ segment number ]
    #Variables are in the order dens, temp, dark , velo - later changing this to a dict would be much more conveniant and elegant



    #Init lists of filaments, and dict to help address this
    vardict = {'dens':0,'kine':1}

    profiles = [ [ [] for x in range(filament_num) ] for i in range(2)]
    
     

    #Determine total list of files
    filelist = sorted(os.listdir(directory))
    #For each file in this list, load data, and append to correct part of profiles structure

   
    filelen = len(filelist)

    denslist = filelist[0:filelen/2]
    kinslist = filelist[filelen/2:]

    stordict = {}
    
    for stor,(i,file_in_dir) in yt.parallel_objects(enumerate(filelist),storage=stordict):
        #prof = yt.load(''.join([directory,'/',file_in_dir]))
        key_id = (vardict[file_in_dir[0:4]],int(file_in_dir[7:10]),int(file_in_dir[13:16]))

        stor.result = yt.load(''.join([directory,'/',file_in_dir]))
        stor.result_id = key_id
        print key_id

    for (key0,key1,key2),value in sorted(stordict):
        profiles[key0][key1].append(values) 

  
    
    return profiles

def load_prof(directory,filament_num,var,keep_list):
    #~~~~~~~~~~~~~~~~#
    # NOT  TESTED    #
    #                #
    #~~~~~~~~~~~~~~~~#

#Function which the profile data for one variable from disk
    #Only returns the array of variable bins, and no other profile metadata - for that use load_all_profs
    #Will also filter out any segments of filaments deemed to be within a halo
    
    #Init list to return

    var_profs = [[] for i in range(filament_num)]

    

    #Gather list of filenames, ergo gathering list of fils and segs
    import os
 
    filelist = sorted(os.listdir(directory))
    #At this stage we can split the filelist into two depending on variable required
    #This is due to density and temp/velocity profiles being stored seperately.
    vardict = {'density':0,'dark_matter_density':0,"temperature":1,"cylindrical_radial_velocity":1}
    profnum = len(filelist) / 2
    filelist = filelist[profnum * vardict[var]:profnum * (vardict[var]+ 1)]
    
    #Parallize the import of the data to speed this process up
    #Use yt's easy to use parallel_objects implementation
    
    storage = {}

    #Gather x data from one profile, incase this is also needed
    x = yt.load(''.join([directory,'/',filelist[0]])).data['x']

    filelist=filelist[:10]
    for sto, file_in_dir in yt.parallel_objects(filelist, storage=storage):
        #Determine filament and segment number
        
        filnum = int(file_in_dir[7:10])
        segnum = int(file_in_dir[13:16])
        
        #Check to see if segment is within a halo, if so, disregard
        #if keep_list[filnum][segnum] == True:
             
        prof = yt.load(''.join([directory,'/',file_in_dir])).data[var]
        sto.result = prof
        sto.result_id = "%d_%d" %(filnum,segnum)
        
        
            

    for (fil,seg),prof in sorted(storage.items()):
        var_prof[fil].append(prof)
        
               
    
    return var_profs,x

                                    



if __name__ == "__main__" :
    
    
    #Load Dataset

    dataset = sys.argv[1]

    ds = yt.load(dataset)
    #halos = generate_halos(ds)
    dsname = dataset[dataset.rfind("/",1):]

    halos,radii,pos = generate_halos(ds,dsname)
  

    filaments = filament.gen_filaments(sys.argv[2])
                        
                               
    keep_arr = sep_fils(filaments,radii,pos,ds)
    
    make_save_profiles(filaments,dsname,ds,check=True)
    

    
    if yt.is_root():

        np.save(''.join(['/shome/mackie/data/',dsname,'/filkeep.npy']),keep_arr)
        np.save(''.join(['/shome/mackie/data/',dsname,'/fils.npy']), filaments)

    print("Completed")
