##!/Users/martinmccullagh/anaconda/bin/python
# NOTE: will have to point to version of python and have MDAnalysis library

# USAGE: dist_windows.py [config file]
# CONFIG FILE FORMAT:
# psffile = [psf, prmtop, gro, or pdb file]
# dcdfile = [traj file in format trr, dcd, etc]
# atom_sel_1 = [CHARMM style atom selection for first group of atoms]
# atom_sel_2 = [CHARMM style atom selection for second group of atoms]
# dist_max = [maximum distance]
# dist_min = [minimum distance]
# dist_delta = [distance bin size]


# load libraries
import sys
import os
import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt




#################################################################################################################
##############################################     SUBROUTINEs     ##############################################
#################################################################################################################


# define subroutines
def ParseConfigFile(cfg_file):
	global dist_min,dist_max,dist_delta,output_hist_name,equilib,counts_out,low_end,r_free_out,u_free_out,pmf_vert
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			print "Option:", option, " Value:", value
			# check value
			if option.lower()=='dist_min':
				dist_min = float(value)
			elif option.lower()=='dist_max':
				dist_max = float(value)
			elif option.lower()=='dist_delta':
				dist_delta = float(value)
			elif option.lower()=='hist_out':
				output_hist_name = value
			elif option.lower()=='r_free_out':
				r_free_out = value
			elif option.lower()=='u_free_out':
				u_free_out = value
			elif option.lower()=='equilibration_steps':
				equilib = int(value)
			elif option.lower()=='bin_counts_outfile':
				counts_out = value
			elif option.lower()=='ignore_bins_below':
				low_end = int(value)
			elif option.lower()=='wham_vertical_correction':
				pmf_vert = value
			else :
				print "Option:", option, " is not recognized"
	
	f.close()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False



#################################################################################################################
##############################################     MAIN PROGRAM     #############################################
#################################################################################################################

ParseConfigFile(sys.argv[1])

k = 0.001987 # kCal K^-1 mol^-1
T = 298.0 # K
kT = k*T
wham_spring_constant = 20.0 # wham spring constant is 2x the AMBER spring constant
restraint_const = 0.5*wham_spring_constant/(kT) ## 1/2 * wham_spring_constant * 1/(kT) [in kcal/mol]
total_bins = int((dist_max - dist_min)/dist_delta) # distance of bin is bin# * dist_delta + dist_min

edges = np.zeros(total_bins+1)
for j in range(total_bins+1):
	edges[j] = dist_min + j*dist_delta

half_bins = np.zeros(total_bins)
for k in range(total_bins):
	half_bins[k] = edges[k] + dist_delta/2.0 # Mid point of each bin.


count = 0
window_value = []
for dirs, subdirs, files in os.walk('./'):
	if dirs[2:4] == "us":
		window = float(dirs[4:])
		window_string = dirs[4:]
		if count == 0 or window < window_min:
			window_min = window
		if count == 0 or window > window_max:
			window_max = window
		check_file = dirs + "/PDI.us." + window_string + ".run01.dat"
		print "Checking " + check_file
		if os.path.isfile(check_file):
			# open file and create histogram
			datalist = np.loadtxt(check_file)


			### Start binning algorithm ###

			nColumns = datalist.shape[1]
			nSteps = datalist.shape[0] # rows
			post_equil_steps = nSteps - equilib # This way I could subtract the 5 ns of equilibration here.
		#	print 'Number of steps: %d, number of steps after equilibration: %d' %(nSteps,post_equil_steps) # convert steps to time so its easier to understand how much time is being taken out 

			weighting_constant = nSteps*dist_delta*4.0*np.pi
			restraint_position = window/10.0 # window put into angstroms
		#	print 'weighting_constant, %f' %weighting_constant

		#	x_min = np.ndarray.min(datalist[equilib:,1]) # the minimum collective variable value.
		#	x_max = np.ndarray.max(datalist[equilib:,1]) # the maximum collective variable value.
		#	bins = int((x_max - x_min)/dist_delta) # divide distance the dat file spans into bins of dist_delta size.
			
		#	print 'Window %d has %d bins' %(window,bins)

			weighted_counts = np.zeros(total_bins)
			counts = np.zeros(total_bins) # will be the number of points in each bin.
			for j in range(equilib,nSteps): # loop inclusive of equilib, exclusive of nSteps.
				index = int((datalist[j][1] - dist_min)/dist_delta) # index data as nearest int steps from dist_min. This makes every window's histogram span the entire range and be positioned correctly wrt other windows
				if index not in range(total_bins):
					print 'Timestep %d is trying to be binned in index %d...' %(j,index)
				
				if index == total_bins: # Data points binned into the last bin don't have a bin to go to so we move them to the last actual bin. Index of -1 is last index. -2 would be second to last..
					weighted_counts[-1] += np.exp(restraint_const * (datalist[j][1] - restraint_position)**2)/(weighting_constant * datalist[j][1]**2) # nSteps*dist_delta ---> normalizes the probability density; 4.0*np.pi*datalist[j][1]**2 ---> volume correction
					counts[-1] += 1
				else:
					weighted_counts[index] += np.exp(restraint_const * (datalist[j][1] - restraint_position)**2)/(weighting_constant * datalist[j][1]**2) # nSteps*dist_delta ---> normalizes the probability density; 4.0*np.pi*datalist[j][1]**2 ---> volume correction
					counts[index] += 1

			for j in range(total_bins):
				if counts[j] <= low_end: # Ignore bins with few points
					weighted_counts[j] = 0

			############## Calculate Restricted Free Energy ###############
			#r_free_energy = np.zeros(total_bins)
			#for j in range(total_bins):
			#	r_free_energy[j] = -kT*np.log(weighted_counts[j])
			#min_fe = np.ndarray.min(r_free_energy)
			#r_free_energy -= min_fe # Make the FE have its minimum at zero

			############## Calculate Unrestricted Free Energy ###############

			u_free_energy = np.zeros(total_bins)
			restraint_position = window/10.0 # window put into angstroms
			#print 'restraint position for window # %d: %f' %(window,restraint_position)
			for j in range(total_bins): # first two u_free_energy calculations are okay if restraint isn't weighted in weighted_counts
				#u_free_energy[j] = r_free_energy[j] - kT * (restraint_const * (half_bins[j] - restraint_position)**2) 
				#u_free_energy[j] = -kT * np.log(weighted_counts[j] * np.exp(restraint_const * (half_bins[j] - restraint_position)**2))
				u_free_energy[j] = -kT * np.log(weighted_counts[j])
			

			##############################
			if count == 0:
				temp_hist = weighted_counts
				bin_counts = counts
				#temp_r_free = r_free_energy
				temp_u_free = u_free_energy

			else:
				temp_hist = np.column_stack((temp_hist,weighted_counts))
				bin_counts = np.column_stack((bin_counts,counts))
				#temp_r_free = np.column_stack((temp_r_free,r_free_energy))
				temp_u_free = np.column_stack((temp_u_free,u_free_energy))

			##############################

			count += 1
			window_value.append(window)


col_var = np.zeros(total_bins)
for t in range(total_bins):
	col_var[t] = dist_min + t*dist_delta
end_hist = np.column_stack((col_var,temp_hist))
np.savetxt(output_hist_name,end_hist)

end_counts = np.column_stack((col_var,bin_counts))
np.savetxt(counts_out,end_counts)

#end_r_free = np.column_stack((col_var,temp_r_free))
#np.savetxt(r_free_out,end_r_free)

pmf_vert_correction = np.loadtxt(pmf_vert)
print temp_u_free.shape
for j in range(len(temp_u_free[1])):
	for i in range(len(temp_u_free[0])):
		temp_u_free[i][j] += pmf_vert_correction[j]

end_u_free = np.column_stack((col_var,temp_u_free))
np.savetxt(u_free_out,end_u_free)


