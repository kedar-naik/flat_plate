#!/usr/bin/env python

# Use this script for validating flat-plate cases from SU^2
# Author: Kedar R. Naik
# Date: June 2014

################################# user inputs #################################

# filename
results_filename = "surface_flow350.dat"
#results_filename = "surface_flow_303.dat"
#results_filename = "flow350.dat"
#results_filename = "flow_303.dat"

# variables of interest (flow.dat and surface_flow.dat)
ref1_variable = 'x'				# search along ref1
ref2_variable = 'y'				# search along ref2

# return corresponding value
val1_variable = "Conservative_1"
val2_variable = "Conservative_2"

val3_variable = "heat_flux"
#val3_variable = "skin_friction_coefficient"	
#val3_variable = "temperature"

# reference1 bounds
ref1_min = 0.97000
ref1_max = 0.97010

#ref1_min = 1.89
#ref1_max = 1.92

# all of the coarse mesh
ref1_min = 0
ref1_max = 0.3

# end of the coarse mesh
#ref1_min = 0.27
#ref1_max = 0.29

#ref1_min = 1.95
#ref1_max = 1.96

# ADD AN INTERPOLATION CAPABILITY?

# reference2 bounds
ref2_min = 0
ref2_max = 4.0123e-3

#ref2_min = 0
#ref2_max = 0.04

# second row from the bottom of coarse mesh
ref2_min = 1.6e-5
ref2_max = 1.7e-5

ref2_min = -1
ref2_max = 1

# axis labels
abscissa_label = 'Re_x'
ordinate_label = 'C_f'

#abscissa_label = r'\theta = \frac{T-T_w}{T_\infty-T_w}'
#ordinate_label = r'\eta'
#
#abscissa_label = r'1-\frac{\rho}{\rho_\infty}'
#ordinate_label = r'\eta'

#abscissa_label = r'\frac{u}{U_\infty}'
#ordinate_label = r'y'

abscissa_label = 'Re_x'
ordinate_label = 'Nu_x'

# data label
#data_label = 'SU2 (137x97)'
data_label = 'SU2 (65x65)'
data_label = r'SU2 ($T_w$ = 350 K)'

# title
#title = 'flat plate, turbulent, SA, coarse mesh: 137x97'
#title = 'flat plate, turbulent, SA, fine mesh: 545x385'
#title = 'flat plate, turbulent, SA, coarse mesh: 137x97, x = 0.970 m'
#title = 'flat plate, turbulent, SA, coarse mesh: 137x97, x = 1.903 m'
#title = 'flat plate, turbulent, SA, fine mesh: 545x385, x = 0.970 m'
#title = 'flat plate, turbulent, SA, fine mesh: 545x385, x = 1.903 m'
#title = 'flat plate, laminar'
title = ''

# plotting axes
custom_axes = 'no'			# ('yes','no')

x_min = -0.02   # for density profiles
x_max = 0.16
y_min = 0
y_max = 8

#x_min = 0       # for temp profiles
#x_max = 1.2
#y_min = 0
#y_max = 8

# Are there points in a two-column file that you wish to plot too?
points_from_file = 'yes'            # ('yes', 'no')

points_filename = 'nu_303.txt'
header_lines = 1
points_label = r'SU2 ($T_w$ = 303 K)'
points_marker = 'b.-'

# Are you computing temperature profiles parallel to the plate?
temp_profiles = 'no'           # ('yes', 'no')
T_w = 350			        # isothermal wall temperature, [K]

# Are you computing density profiles parallel to the plate?
density_profiles = 'no'        # ('yes', 'no')

# Would you like to save the figure? What should we call it?
save_plot = 'yes'			# ('yes','no')
save_pic_as = 'current_plot.png'

# Have additional data sets (e.g. experimental data) been defined below?
more_data = 'yes'			# ('yes','no')

# Are you going to alter the data coming from the output file yourself below?
altered_data = 'yes'			# ('yes','no')

# Would you like a file to be rewritten containing the extacted/plotted data?
write_file = 'yes'			# ('yes','no')
save_file_as = 'extracted_points.txt'

###############################################################################
import matplotlib.pyplot as plt
import math
#import numpy as np

plt.close("all")

# plot the SU^2 results

# open the results file
results_file = open(results_filename,'r')

# search along column
ref1_variable = ref1_variable.lower()
ref2_variable = ref2_variable.lower()
val1_variable = val1_variable.lower()
val2_variable = val2_variable.lower()
val3_variable = val3_variable.lower()

# go through the header appropriately
line_counter = 1
val_counter = 0
reference1 = []
reference2 = []
values1 = []
values2 = []
values3 = []
for line in results_file:
  
  if line_counter == 2:
    
    # read the variable names
    stripped_line = line.strip('VARIABLES = "').lower()
    stripped_line = stripped_line.lower()
    stripped_line = stripped_line.rstrip('"\n')
    stripped_line = "".join(stripped_line.split('"'))
    stripped_line = "".join(stripped_line.split())
    variables = stripped_line.split(',')
    
    # figure out which column contains the reference value
    entry_counter = 1
    for entry in variables:
      if entry == ref1_variable:
        ref1_column = entry_counter
      if entry == ref2_variable:
        ref2_column = entry_counter
      if entry == val1_variable:
        val1_column = entry_counter
      if entry == val2_variable:
        val2_column = entry_counter
      if entry == val3_variable:
        val3_column = entry_counter
      entry_counter += 1
    
    # print the column numbers to the screen
    print "variable ", ref1_variable, 'is found in column ', ref1_column
    print "variable ", ref2_variable, 'is found in column ', ref2_column
    print "variable ", val1_variable, 'is found in column ', val1_column
    print "variable ", val2_variable, 'is found in column ', val2_column
    print "variable ", val3_variable, 'is found in column ', val3_column

    # convert between column no. and index
    ref1_index = ref1_column-1
    ref2_index = ref2_column-1
    val1_index = val1_column-1
    val2_index = val2_column-1
    val3_index = val3_column-1

  elif line_counter == 3:
    
    # read and print the total number of nodes
    line = line.split()
    nodes = int(line[2].rstrip(','))
    print "nodes in file = ", nodes
 
  elif line_counter >= 4:
    
    if line_counter >= 4 and line_counter <= nodes+3:     
 
      # tokenize the line
      line = line.split()
      
      # recast, rename
      line_at_ref1 = float(line[ref1_index])
      line_at_ref2 = float(line[ref2_index])
      line_at_val1 = float(line[val1_index])
      line_at_val2 = float(line[val2_index])
      line_at_val3 = float(line[val3_index])
      
      # check range of the first refernce value
      if line_at_ref1 >= ref1_min and line_at_ref1 <= ref1_max:
         
        # check range of the second reference value
        if line_at_ref2 >= ref2_min and line_at_ref2 <= ref2_max:
          
          # record the quantity of interest
          reference1.append(line_at_ref1)
          reference2.append(line_at_ref2)
          values1.append(line_at_val1)
          values2.append(line_at_val2)
          values3.append(line_at_val3)
          val_counter += 1
          
      #print values1

  else:
    pass
  
  line_counter += 1
  
print "number of values recorded = ", val_counter  

  
###############################################################################
def two_col_read(filename, header_rows):
  """
  This function reads two columns of data from a file, ignores the header, 
  and returns the two columns as lists.
  Input: The name of the file to be read, the number of rows in the header
  Output: The two columns as lists.
  """
  our_file = open(filename,'r')
  lines = our_file.readlines()
  x = []
  y = []
  line_counter = 1
  for line in lines:
    if line_counter > header_rows:
      split_line = line.split()
      if len(split_line) != 0:
        x.append(float(split_line[0]))
        y.append(float(split_line[1]))
    line_counter += 1
  return (x,y)
  
###############################################################################
def pohlhausen_profile(Pr):
    """
    Given the Prandtl number of the flat-plate flow, this function returns the
    nondimensional temperature profile of Pohlhausen (1921), theta, along with 
    the corresponding values of nondimensional wall distance, eta.
    N.B. theta = (T-T_w)/(T_inf - T_w)
    Input: Pr
    Output: theta and eta as lists
    """
    
    from scipy.integrate import simps
    
    # eta and Blasius solution values from Rochester lecture notes
    # (N.B. The values given in White are different, possibly wrong)
    eta = [0 , 0.25, 0.50, 0.75, 1.00,
               1.25, 1.50, 1.75, 2.00,
               2.25, 2.50, 2.75, 3.00,
               3.25, 3.50, 3.75, 4.00,
               4.25, 4.50, 4.75, 5.00,
               5.25, 5.50, 5.75, 6.00,
               6.25, 6.50, 6.75, 7.00,
               7.25, 7.50, 7.75, 8.00]
               
    f = [0, 0.0104, 0.0415, 0.0933, 0.1656,
             0.2580, 0.3701, 0.5011, 0.6500,
             0.8156, 0.9963, 1.1906, 1.3968,
             1.6131, 1.8377, 2.0691, 2.3057,
             2.5464, 2.7901, 3.0360, 3.2833,
             3.5316, 3.7806, 4.0300, 4.2796,
             4.5294, 4.7793, 5.0293, 5.2792,
             5.5292, 5.7792, 6.0292, 6.2792]

    dfdn = [0, 0.0830, 0.1659, 0.2483, 0.3298,
                0.4096, 0.4868, 0.5605, 0.6298,
                0.6936, 0.7513, 0.8022, 0.8460,
                0.8829, 0.9130, 0.9370, 0.9555,
                0.9694, 0.9795, 0.9867, 0.9915,
                0.9948, 0.9969, 0.9982, 0.9990,
                0.9994, 0.9997, 0.9998, 0.9999,
                1.0000, 1.0000, 1.0000, 1.0000]
                
    d2fdn2 = [0.3321, 0.3319, 0.3309, 0.3282, 0.3230,
                      0.3146, 0.3026, 0.2866, 0.2668,
                      0.2434, 0.2174, 0.1897, 0.1614,
                      0.1337, 0.1078, 0.0844, 0.0642,
                      0.0474, 0.0340, 0.0236, 0.0159,
                      0.0104, 0.0066, 0.0040, 0.0024,
                      0.0014, 0.0008, 0.0004, 0.0002,
                      0.0001, 0.0001, 0.0000, 0.0000]
    
    # take the integrals given in the MIT lecture notes
    N_values = len(eta)
    theta = [float('nan') for x in range(N_values)]
    for i in range(N_values):
        
        integrand = [pow(entry,Pr) for entry in d2fdn2[i:N_values]]
        numerator = simps(integrand,eta[i:N_values])
        
        integrand = [pow(entry,Pr) for entry in d2fdn2[0:N_values]]
        denominator = simps(integrand,eta[0:N_values])
        
        theta[i] = 1 - numerator/denominator
    
    return (eta,theta)

######################## alter references and/or values #######################

# hardcoded values
a = 347.202 				# sound speed, [m/s]
Ma_inf = 0.1				# freestream mach number
mu = 1.84492e-5				# dynamic viscosity, [N.s/m^2]
U_inf = 69.5429				# freestream velocity, [m/s]
U_inf = 34.7715
T_inf = 300				# freestream temperature, [K]
rho_inf = 1.13753			# freestream density, [kg/m^3]
#rho_inf = 2.27506
cp = 1005				# specific heat capacity, [J/kg.K]
Pr = 0.72				# Prandtl number
k = 0.02618                 # thermal conductivity, [W/m.K]

# gotten from top-right corner of volume flow
rho_inf = 1.15419   # T_w = 350
rho_inf = 1.155    # [kg/m^3]T_w = 303
#T_inf = 299.4012    # [K] T_w = 303 and T_w = 350

# Let T_inf be defined by the temperature found at the of the outlet boundary
# (N.B. This assumes that values3 = temperatures along the outlet, starting 
#  from the top fo the mesh and going down)
if temp_profiles == 'yes':
    T_inf = values3[0]

# similarly for density
if density_profiles == 'yes':
    rho_inf = values1[0]
    
# derived quantities
#U_inf = Ma_inf*a		
r = pow(Pr,1/2)			# recovery factor, laminar flow
#r = 0.84771                 # white, pg. 515
T_aw = T_inf + r*pow(U_inf,2)/(2*cp)	# adibatic-wall temperature, [K]

# variable declaration
u_over_U = []
Re_x = []
Cf_Blasius = []
nondim_T = []
eta = []
nondim_T_CB = []
nondim_T2 = []
eta_over_2 = []
Nu_x_Blasius = []
Nu_x = []

# renaming reported quantities
for value in range(val_counter):
  
  rho = values1[value]
  rhou = values2[value]
  u = rhou/rho
  u_over_U.append(u/U_inf)
  
  # rename the references
  x = reference1[value]
  y = reference2[value]

  # rename the values
  T = values3[value]
  q_w = values3[value]        # heat flux at the wall
  #rho = values2[value]
  #rhou = values3[value]
  #u = rhou/rho
  #u_over_U.append(u/U_inf)
  
  
  # perform computations

  # crocco-busemann relation (laminar, compressible bounardary layers)
  # N.B. For this to work properly, you need to be computing u 
  # Need: rho=values1, and rhou = value2 
  T_CB = T_w + (T_aw-T_w)*(u/(U_inf)) - r*pow(u,2)/(2*cp)
  nondim_T_CB.append((T_CB - T_w)/(T_inf - T_w)) 
  
  # local Nusselt number from simulation  
  Nu_x.append(-q_w*x/(k*(T_w-T_inf)))
  
  # local Reynolds number
  nu = mu/rho
  nu = mu/rho_inf
  Re_x.append(U_inf*x/nu)
  
  Cf_Blasius.append(0.6641/math.sqrt(Re_x[value]))
  
  nondim_T.append((T-T_w)/(T_inf-T_w))
  eta.append(y/math.sqrt(nu*x/U_inf))
  
  nondim_T2.append((T-T_inf)/(T_w-T_inf))
  eta_over_2.append((y/2)/math.sqrt(nu*x/U_inf))
  
  eta_Pohl, nondim_T_Pohl = pohlhausen_profile(Pr)
  
  Nu_x_Blasius.append(0.332*pow(Re_x[value],0.5)*pow(Pr,float(1.0/3.0)))
  
# set the values that are to be plotted
abscissas = u_over_U 
ordinates = reference2

# for skin friction
abscissas = Re_x
ordinates = values3

#abscissas = nondim_T
#ordinates = eta
#
#x = reference1
#y = reference2
#rho = values1
#one_minus_nondim_rho = [1-(entry/rho_inf) for entry in rho]
#abscissas = one_minus_nondim_rho
#ordinates = eta

#abscissas = eta_over_2
#ordinates = nondim_T2

abscissas = Re_x
ordinates = Nu_x

###############################################################################
###############################################################################
# add more data sets
# plot your additional data sets or altered values as you like

# skin friction
#x_cfl3d, Cf_cfl3d = two_col_read("cfl3d.dat",2)
#x_fun3d, Cf_fun3d = two_col_read("fun3d.dat",2)
#plt.plot(x_cfl3d,Cf_cfl3d,'b-',label="CFL3D")
#plt.plot(x_fun3d,Cf_fun3d,'r-',label="FUN3D")

# normalized velocity in x
#u_cfl3d_097, y_cfl3d_097 = two_col_read("y_vs_u_cfl3d_0.970.dat",2)
#u_cfl3d_190, y_cfl3d_190 = two_col_read("y_vs_u_cfl3d_1.903.dat",2)
#plt.plot(u_cfl3d_097,y_cfl3d_097,'b-',label="CFL3D")
#plt.plot(u_cfl3d_190,y_cfl3d_190,'b-',label="CFL3D")

# the blasius solution for skin friction
#plt.plot(abscissas, Cf_Blasius, 'k-', label="Blasius")

# the blasius solution for local Nusselt number
plt.plot(abscissas, Nu_x_Blasius, 'k-', label="Blasius")

# plot the pohlhausen solution
#plt.plot(nondim_T_Pohl, eta_Pohl,'k-',label="Pohlhausen")

# plot crocco-busemann profile
#plt.hold(True)
#plt.plot(nondim_T_CB, eta,'g-',label="Crocco-Busemann")


plt.hold(True)

##############################################################################
# plot the contents of a two-column file
if points_from_file == 'yes':
    points_x, points_y = two_col_read(points_filename, header_lines)
    plt.plot(points_x, points_y,points_marker,label=points_label)
    plt.hold(True)
###############################################################################

# generate plot
plt.rc('text', usetex=True)			# for using latex
plt.rc('font',family='serif')			# setting font

if altered_data == 'no':
  abscissas = reference1
  ordinates = values1
  if more_data == 'no':
    plt.hold(False)
plt.plot(abscissas,ordinates,'r.-',label=data_label)
plt.xlabel('$'+abscissa_label+'$',fontsize=16)
plt.ylabel('$'+ordinate_label+'$',fontsize=16)

# finishing
plt.grid(True)
plt.title(title)
if custom_axes == 'yes':
  plt.axis((x_min,x_max,y_min,y_max))
plt.legend(loc = 'upper center')
if more_data == 'yes':
  plt.hold(False)
fig1 = plt.gcf()
plt.show()

# save figure
if save_plot == 'yes':
  fig1.savefig(save_pic_as, dpi=200)
  print 'The plot has been saved as:',save_pic_as
 
# generate file of extracted/plotted results
if write_file == 'yes':
  file = open(save_file_as,'w')
  file.write(abscissa_label+'\t'+ordinate_label+'\n')
  for value in range(val_counter):
    file.write(str(abscissas[value])+'\t'+str(ordinates[value])+'\n')
  file.close()
  print 'The extracted profiles have been saved to a file called:',save_file_as

