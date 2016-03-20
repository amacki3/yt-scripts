import numpy as np
import yt
import yt.units as u

def get_cyl_params(start,end):
	#Function to determine the normal,height and midpoint for the cylinder object given a start and end point
	start = start[0:3]
	end = end[0:3]
	#Determine midpoint of filament and orientation and normal of cylinder
	norm = []
	midpoint = []
	for s,e in zip(start,end):
		if abs(s-e) > 0.5:
			norm.append( 1 + e -s)
			midpoint.append( (s+e)/2.0 - 0.5)
		else:
			norm.append(s-e)
			midpoint.append((s+e)/2.0)
	normal = np.array(norm)
	height = np.sqrt(normal.dot(normal)) / 2.0
	del norm
	return normal,height,midpoint

def gen_cyl(start,end,radius,ds,bulk_vel=True):
	#Function to create a cylinder object given a start point,end point and radius.
	#This also ensures the correct bulk velocity, if required, is set for this cyliner object
	normal,height,midpoint = get_cyl_params(start,end)
	
        #Create Cylinder Object along filament, radius equal to radius passed in MPc
	try:
		radius = radius.in_units("Mpc")
	except:
		radius = radius * u.Mpc

	
	
	cyl = ds.disk(midpoint,normal,radius=(radius),height=height)
	
	if bulk_vel == True:
        #Set bulk velocity for cylinder
		cyl2 = ds.disk(midpoint,normal,radius=(radius),height=height)
		bulk_vel = cyl2.quantities.bulk_velocity()
		del cyl2
		cyl.set_field_parameter("bulk_velocity", bulk_vel)
	
	return cyl

def det_fil_profile(start,end,var,weight,radius,ds):
	
	#Function to return a list of the weighted average of a variable for a list of concentric disks

	#Get the disk object

	cyl = gen_cyl(start,end,radius,ds)

	#Create a profile
	if weight == "":
		weight = None
	profile = yt.create_profile(cyl,"cylindrical_radius", var, weight_field=weight,n_bins=10,units={"cylindrical_radius":'Mpc'},extrema={"cylindrical_radius":(0.1,2.0)})
	
	return profile



	
	


