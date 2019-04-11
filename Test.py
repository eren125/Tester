#######################################################################
############Treating DL_Poly energies, forces, coords data.############
#######################################################################

import numpy as np

n = 4
natm = 6
List = open("out.txt").readlines()  ## Creating a list with the output file

def energies(n,natm):
 	E = [ [0 for k in range(natm)] for l in range(n)]
# adding energies in the energy table; at the end we will have the time dependence and atom dependence in one table.
	s = 0
	for i in range(n): # i is time
		s += 1 # skipping the first line and line with total energy
		for j in range(natm): # j is atom number
			# we split the string into a list of stings and the velue of energy is contained in the 4th. We change it into a float, so we can do further calculations.
			E[i][j] = float( List[s].split()[4] ) 
			s += 3 # skipping coords and forces lines
	return E

#######################################################################################################################################################################################

def forcesx(n,natm):
 	E = [ [0 for k in range(natm)] for l in range(n)]
# adding energies in the energy table; at the end we will have the time dependence and atom dependence in one table.
	s = 1
	for i in range(n): # i is time
		s += 1 # skipping the first line and line with total energy
		for j in range(natm): # j is atom number
			# we split the string into a list of stings and the velue of energy is contained in the 4th. We change it into a float, so we can do further calculations.
			E[i][j] = float( List[s].split()[4] ) 
			s += 3 # skipping coords and forces lines
	return E

#######################################################################################################################################################################################

def forcesy(n,natm):
	E = [ [0 for k in range(natm)] for l in range(n)]
# adding energies in the energy table; at the end we will have the time dependence and atom dependence in one table.
	s = 1
	for i in range(n): # i is time
		s += 1 # skipping the first line and line with total energy
		for j in range(natm): # j is atom number
			# we split the string into a list of stings and the velue of energy is contained in the 4th. We change it into a float, so we can do further calculations.
			E[i][j] = float( List[s].split()[5] ) 
			s += 3 # skipping coords and forces lines
	return E

#######################################################################################################################################################################################

def forcesz(n,natm):
	E = [ [0 for k in range(natm)] for l in range(n)]
# adding energies in the energy table; at the end we will have the time dependence and atom dependence in one table.
	s = 1
	for i in range(n): # i is time
		s += 1 # skipping the first line and line with total energy
		for j in range(natm): # j is atom number
			# we split the string into a list of stings and the velue of energy is contained in the 4th. We change it into a float, so we can do further calculations.
			E[i][j] = float( List[s].split()[6] ) 
			s += 3 # skipping coords and forces lines
	return E

#######################################################################################################################################################################################

def coordsx(n,natm):
	E = [ [0 for k in range(natm)] for l in range(n)]
# adding energies in the energy table; at the end we will have the time dependence and atom dependence in one table.
	s = 2
	for i in range(n): # i is time
		s += 1 # skipping the first line and line with total energy
		for j in range(natm): # j is atom number
			# we split the string into a list of stings and the velue of energy is contained in the 4th. We change it into a float, so we can do further calculations.
			E[i][j] = float( List[s].split()[4] ) 
			s += 3 # skipping coords and forces lines
	return E

#######################################################################################################################################################################################

def coordsy(n,natm):
	E = [ [0 for k in range(natm)] for l in range(n)]
# adding energies in the energy table; at the end we will have the time dependence and atom dependence in one table.
	s = 2
	for i in range(n): # i is time
		s += 1 # skipping the first line and line with total energy
		for j in range(natm): # j is atom number
			# we split the string into a list of stings and the velue of energy is contained in the 4th. We change it into a float, so we can do further calculations.
			E[i][j] = float( List[s].split()[5] ) 
			s += 3 # skipping coords and forces lines
	return E

#######################################################################################################################################################################################

def coordsz(n,natm):
	E = [ [0 for k in range(natm)] for l in range(n)]
# adding energies in the energy table; at the end we will have the time dependence and atom dependence in one table.
	s = 2
	for i in range(n): # i is time
		s += 1 # skipping the first line and line with total energy
		for j in range(natm): # j is atom number
			# we split the string into a list of stings and the velue of energy is contained in the 4th. We change it into a float, so we can do further calculations.
			E[i][j] = float( List[s].split()[6] ) 
			s += 3 # skipping coords and forces lines
	return E

#######################################################################################################################################################################################

def deltaEnergies(n,natm):
	E = energies(n,natm)
	DE = []
	for j in range(natm): # For the atom j,
		de = [] 
		for i in range(n-1): # we put the difference of energy between to two time steps (t_i , t_i+1)
			de.append(E[i+1][j]-E[i][j])
		DE.append(de)	
	return DE

#######################################################################################################################################################################################

def taylorExpansion(n, natm):
	Fx = forcesx(n,natm)
	Fy = forcesy(n,natm)
	Fz = forcesz(n,natm)
	Cx = coordsx(n,natm)
	Cy = coordsy(n,natm)
	Cz = coordsz(n,natm)
	S = []
	for j in range(natm): # For the atom j,
		s = [] 
		for i in range(n-1): # we put the difference of energy between to two time steps (t_i , t_i+1)
			s.append(-(((Cx[i+1][j]-Cx[i][j])*Fx[i][j])+((Cy[i+1][j]-Cy[i][j])*Fy[i][j])+((Cz[i+1][j]-Cz[i][j])*Fz[i][j])))
		S.append(s)	
	return S 

for i in range(1,n):
	if float( List[(3*natm+1)*i].split()[4] ) < float( List[(3*natm+1)*(i+1)].split()[4] ) :
		print("Forces are in the wrong direction between steps", i, (i+1)) 
		print( float( List[(3*natm+1)*i].split()[4] ), float( List[(3*natm+1)*(i+1)].split()[4] ) )
		print("")

###################################################################
## If you want to print the raw values of energies forces coords.##
###################################################################

print("List of energies : ")
print("")
print(np.array(energies(n,natm)))
print("")
print("List of forces : ")
print("")
print(np.array(forcesx(n,natm)))
print(np.array(forcesy(n,natm)))
print(np.array(forcesz(n,natm)))
print("")
print("List of coords :")
print("")
print(np.array(coordsx(n,natm)))
print(np.array(coordsy(n,natm)))
print(np.array(coordsz(n,natm)))
print("")

###############################################################################
## If you want to compare the difference of energy with its Taylor expansion ##
###############################################################################

print("Delta E :")
print("")
print(np.array(deltaEnergies(n,natm)))
print("")
print("Taylor expansion :")
print("")
print(np.array(taylorExpansion(n,natm)))
print("")

