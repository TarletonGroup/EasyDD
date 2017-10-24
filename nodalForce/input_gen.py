# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 15:16:40 2017

@author: daniel
"""

import numpy as np

n_dln = 5000
n_se = 5000
threads_per_block = 128 # Max 256 for this graphics card

f = open('cuda_input.txt', 'w')
f.write('%d\n' %n_dln) # N dislocations
f.write('%d\n' %n_se) # N surface elements
f.write('%d\n' %threads_per_block) # N surface elements
f.write('\n')
for i in range(n_dln): # Dislocation nodes
    nodes = np.random.random(3)
    f.write('%f %f %f\n' %(nodes[0], nodes[1], nodes[2]))
    f.write('%f %f %f\n' %(nodes[0]+1, nodes[1]+1, nodes[2]+1))
f.write('\n')
for i in range(n_se): # Surface element nodes
    nodes = np.random.random(3)
    f.write('%f %f %f\n' %(nodes[0], nodes[1], nodes[2]))
    f.write('%f %f %f\n' %(nodes[0]+1, nodes[1], nodes[2]))
    f.write('%f %f %f\n' %(nodes[0], nodes[1]+1, nodes[2]))
    f.write('%f %f %f\n' %(nodes[0]+1, nodes[1]+1, nodes[2]))
    #f.write('%f %f %f\n' %(nodes[0], nodes[1], nodes[2]))
f.write('\n')
for i in range(n_dln): # Burgers' vectors
    nodes = np.random.random(3)
    f.write('%f %f %f\n' %(nodes[0], nodes[1], nodes[2]))
f.write('\n')
nodes = np.random.random(3)
f.write('%f\n' %nodes[0]) # mu
f.write('%f\n' %nodes[1]) # nu
f.write('%f\n' %(nodes[2]*0.1)) # nu
f.close()
