import yt
import matplotlib.pyplot as plt
import filament
import numpy as np
import sys
import math

def projplot(ds,filaments,axis):
	


	axes = {'x':0,'y':1,'z':2}
	x = axes[axis]
	

       	plot = yt.ProjectionPlot(ds,axis, "density",weight_field="density")
	
	#Gather line segments
	for my_fil in filaments:
		for i in xrange(1,len(my_fil)):
			y1 = my_fil[i-1][(x+1)%3]
		       	y2 = my_fil[i][(x+1)%3]
	       		z1 = my_fil[i-1][(x+2)%3]
			z2 = my_fil[i][(x+2)%3]
	#Plot each individual segment
			if (abs(y2-y1) < 0.9) and (abs(z2-z1) < 0.9): 
				plot.annotate_line((y1,z1),(y2,z2),coord_system='axis')
			del y1,y2,z1,z2
	
	return plot

if __name__ =="__main__":
	#load dataset from sys arg 1
	dataset = sys.argv[1]
	ds = yt.load(dataset)

	#load filaments
	try:
		filaments = filament.load_filaments(sys.argv[2])
	except IOError:
		filaments = filament.gen_filaments(sys.argv[2])
	print("Filaments read.")
	
	plotted = projplot(ds,filaments,'x')
	plotted.save("".join([sys.argv[3],'x.png']))
	'''
	plotted = projplot(ds,filaments,'y')
        plotted.save("".join([sys.argv[3],'y.png']))
	plotted = projplot(ds,filaments,'z')
        plotted.save("".join([sys.argv[3],'z.png']))
	'''
