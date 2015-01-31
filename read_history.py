# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 18:10:38 2015

@author: Kedar
"""
# filename
results_filename = "history350.dat"

# variables of interest
variable = "Res_Flow[0]"               

# axis labels
abscissa_label = r'iteration'
ordinate_label = r'\log_{10}(R_{\rho})'

# data label
data_label = 'SU2'

# title
title = ''

# Would you like to save the figure? What should we call it?
save_plot = 'no'			# ('yes','no')
save_pic_as = 'history_plot.png'

# Would you like a file to be rewritten containing the extacted/plotted data?
write_file = 'no'			# ('yes','no')
save_file_as = 'history.txt'

###############################################################################
import matplotlib.pyplot as plt

plt.close("all")

# open the results file
results_file = open(results_filename,'r')

# make variable to be miniscules
variable = variable.lower()

# go through the header appropriately
line_counter = 1
val_counter = 0
values = []
for line in results_file:
  
  if line_counter == 2:
    
    # read the variable names
    stripped_line = line.strip('VARIABLES = ')
    stripped_line = stripped_line.lower()
    stripped_line = stripped_line.rstrip('"\n')
    stripped_line = stripped_line.lstrip('"')
    variables = stripped_line.split('","')
    
    # figure out which column contains the reference value
    entry_counter = 1
    for entry in variables:
      if entry == variable:
        val_column = entry_counter
      entry_counter += 1
    
    # print the column numbers to the screen
    print "variable ", variable, 'is found in column ', val_column

    # convert between column no. and index
    val_index = val_column-1
 
  elif line_counter >= 4:
    
    # tokenize the line
    line = line.split(',')
  
    # recast, rename
    line_at_val = float(line[val_index])
    
    # store
    values.append(line_at_val)
    
    # number of values recorded
    val_counter += 1
      
  else:
    pass
  
  line_counter += 1
  
print "number of values recorded = ", val_counter  

# start plotting
abscissas = range(val_counter)
ordinates = values

# generate plot
plt.rc('text', usetex=True)			# for using latex
plt.rc('font',family='serif')			# setting font

plt.plot(abscissas,ordinates,'k-',label=data_label)
plt.xlabel('$'+abscissa_label+'$',fontsize=16)
plt.ylabel('$'+ordinate_label+'$',fontsize=16)

# finishing
plt.grid(True)
plt.title(title)

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
  print 'The extracted history has been saved to a file called:',save_file_as

