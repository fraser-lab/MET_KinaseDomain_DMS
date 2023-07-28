# MET_ensemble_alignment_2.pml
# By G.Estevam
# This is a program that opens MET kinase structures, aligns them to nlobe residues, measures their distance, and retrieves additional strcture info

#### pseudocode:
# fetch all pdbs, align, and make them aesthstically pleasing
# create a dictionary for each residue
# interate state so that the residue coordinates for each state are saved in the dictionary for each structure/ object
# run a loop that access the residues coordinate array, assigns each coordinate a variable, and does math to determine the distance between residues
# grab additional structure details: resolution, ligand, resolved JM (based on defined residues)
# append the structural information into a csv file (PDB, residue distances, resolution, ligand, resolved JM)


# reinitialize pymol
reinit

# import modules
import os
import re
import sys
import subprocess
import math
import csv


# structure retrival and aesthetics
#fetch 1R0P 1R1W 2G15 2RFN 2RFS 2WD1 2WGJ 2WKM 3A4P 3CCN 3CD8 3CE3 3CTH 3CTJ 3DKC 3DKF 3DKG 3EFJ 3EFK 3F66 3F82 3I5N 3L8V 3LQ8 3Q6U 3Q6W 3QTI 3R7O 3RHK 3U6H 3U6I 3VW8 3ZBX 3ZC5 3ZCL 3ZXZ 3ZZE 4AOI 4AP7 4DEG 4DEH 4DEI 4EEV 4GG5 4GG7 4IWD 4KNB 4MXC 4R1V 4R1Y 4XMO 4XYF 5DG5 5EOB 5EYC 5EYD 5HLW 5HNI 5HO6 5HOA 5HOR 5HTI 5T3Q 5UAB 5UAD 5YA5 6SD9 6SDC 6SDD 6SDE
fetch 1R0P 1R1W 2G15 2RFN 2RFS 2WD1 2WGJ 2WKM 3A4P 3CCN 3CD8 3CE3 3CTH 3CTJ 3DKC 3DKF 3DKG 3EFJ 3EFK 3F66 3F82 3I5N 3L8V 3LQ8 3Q6U 3Q6W 3QTI 3R7O 3RHK 3U6H 3U6I 3VW8 3ZBX 3ZC5 3ZCL 3ZXZ 3ZZE 4AOI 4AP7 4DEG 4DEH 4DEI 4EEV 4GG5 4GG7 4IWD 4KNB 4MXC 4R1V 4R1Y 4XMO 4XYF 5DG5 5EOB 5EYC 5EYD 5HLW 5HNI 5HO6 5HOA 5HOR 5HTI 5T3Q 5UAB 5UAD 5YA5 6SD9 6SDC 6SDD 6SDE 6UBW 7B3Z 7B3Q 7B42 7B41 7B44 7B43 7B3T 7B3W 7B3V 7B40 7V3R 7V3S 7Y4T 7Y4U 8GVJ 8AN8 8ANS 8OUU 8OUV 8OU7 8OVZ 8OW3 8OWG
sele model 1R0P, resi 1048-1152
remove not chain A
hide everything
show cartoon
sele nlobe, resi 1048-1152
set_color mermaid,[8,144,153]
#set_color plum,[124,29,111]
color white, all
color orange, resi 1025-1070
#color plum, resi 1222-1245 # color activation loop
color mermaid, resi 1117-1134 # color cHelix
set cartoon_fancy_helices, 1
bg_color white
set ray_shadows,0


# align structures
# the goal here is to align all structres to the nlobe sequence
# here the program is aligning each structure individually to the specified nlobe sequence
# things that are important to consider are that we only want to align the same chain (i.e chain A)- some structures contain multiple chains (eventually I want to align these too)
# the approach to do this here is by running a loop that porvides inputs on the command line

python
#pdbs = ['1R0P', '1R1W', '1G15', '2RFN', '2RFS', '2WD1', '2WGJ', '2WKM', '3A4P', '3CCN', '3CD8', '3CE3', '3CTH', '3CTJ', '3DKC', '3DKF', '3DKG', '3EFJ', '3EFK', '3F66', '3F82', '3I5N', '3L8V', '3LQ8', '3Q6U', '3Q6W', '3QTI', '3R7O', '3RHK', '3U6H', '3U6I', '3VW8', '3ZBX', '3ZC5', '3ZCL', '3ZXZ', '3ZZE', '4AOI', '4AP7', '4DEG', '4DEH', '4DEI', '4EEV', '4GG5', '4GG7', '4IWD', '4KNB', '4MXC', '4R1V', '4R1Y', '4XMO', '4XYF', '5DG5', '5EOB', '5EYC', '5EYD', '5HLW', '5HNI', '5HO6', '5HOA', '5HOR', '5HTI', '5T3Q', '5UAB', '5UAD', '5UAF', '5YA5', '6SD9', '6SDC', '6SDD', '6SDE]

pdbs = ['1R0P', '2G15', '1R1W', '2RFN', '2RFS', '2WD1', '2WGJ', '2WKM','3A4P', '3CCN', '3CD8','3CE3', '3CTH', '3CTJ', '3DKC', '3DKF', '3DKG', '3EFJ', '3EFK', '3F66', '3F82', '3I5N','3L8V', '3LQ8', '3Q6U', '3Q6W', '3QTI', '3R7O', '3RHK', '3U6H', '3U6I', '3VW8', '3ZBX', '3ZC5', '3ZCL', '3ZXZ', '3ZZE', '4AOI', '4AP7', '4DEG', '4DEH', '4DEI', '4EEV', '4GG5', '4GG7', '4IWD', '4KNB', '4MXC', '4R1V', '4R1Y', '4XMO', '4XYF', '5DG5', '5EOB', '5EYC', '5EYD', '5HLW', '5HNI', '5HO6', '5HOA', '5HOR', '5HTI', '5T3Q', '5UAB', '5UAD', '5UAF', '5YA5', '6SD9', '6SDC', '6SDD', '6SDE', '6UBW', '7B3Z', '7B3Q', '7B42', '7B41', '7B44', '7B43', '7B3T', '7B3W', '7B3V', '7B40', '7V3R', '7V3S', '7Y4T', '7Y4U', '8AN8','8ANS', '8OUU', '8OUV', '8OU7', '8OVZ', '8OW3', '8OWG']

j = 0
while j<100:
  #cmd.remove()
  #cmd.select('nlobe','resi 1048-1152')
  cmd.align(pdbs[j],'model 1R0P, resi 1048-1152')
  #cmd.align('4IWD','nlobe')
  print (j)
  j = j+1

python end

#pdbs = ['1R0P', '1R1W', '1G15', '2RFN', '2RFS', '2WD1', '2WGJ', '2WKM', '3A4P', '3CCN', '3CD8', '3CE3', '3CTH', '3CTJ', '3DKC', '3DKF', '3DKG', '3EFJ', '3EFK', '3F66', '3F82', '3I5N', '3L8V', '3LQ8', '3Q6U', '3Q6W', '3QTI', '3R7O', '3RHK', '3U6H', '3U6I', '3VW8', '3ZBX', '3ZC5', '3ZCL', '3ZXZ', '3ZZE', '4AOI', '4AP7', '4DEG', '4DEH', '4DEI', '4EEV', '4GG5', '4GG7', '4IWD', '4KNB', '4MXC', '4R1V', '4R1Y', '4XMO', '4XYF', '5DG5', '5EOB', '5EYC', '5EYD', '5HLW', '5HNI', '5HO6', '5HOA', '5HOR', '5HTI', '5T3Q', '5UAB', '5UAD', '5UAF', '5YA5', '6SD9', '6SDC', '6SDD', '6SDE]
#pdbs = ['3Q6W', '3R7O', '4IWD']

resi_1066_xyz_dictionary={}                                                     # creates a dictionary for residue coordinates
resi_1129_xyz_dictionary={}
resi_1062_xyz_dictionary={}
resi_1125_xyz_dictionary={}
resi_1058_xyz_dictionary={}
resi_1121_xyz_dictionary={}

iterate_state 1, (resi 1066 and n. CA),resi_1066_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for V1069 for all objects
iterate_state 1, (resi 1129 and n. CA),resi_1129_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for K1131 for all objects
iterate_state 1, (resi 1062 and n. CA),resi_1062_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for V1069 for all objects
iterate_state 1, (resi 1125 and n. CA),resi_1125_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for K1131 for all objects
iterate_state 1, (resi 1058 and n. CA),resi_1058_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for V1069 for all objects
iterate_state 1, (resi 1121 and n. CA),resi_1121_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for K1131 for all objects


resi_1124_xyz_dictionary={}                                                     # creates a dictionary for residue coordinates
resi_1153_xyz_dictionary={}
resi_1112_xyz_dictionary={}
iterate_state 1, (resi 1124 and n. CA),resi_1124_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for V1069 for all objects
iterate_state 1, (resi 1153 and n. CA),resi_1153_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for K1131 for all objects
iterate_state 1, (resi 1112 and n. CA),resi_1112_xyz_dictionary[model]=[x,y,z]  # appends the coordinates to the dictionary for V1069 for all objects
#print resi_1069_xyz_dictionary


# calculate the residue-to-residue distance for each MET structure

python

pdbs = ['1R0P', '2G15', '1R1W', '2RFN', '2RFS', '2WD1', '2WGJ', '2WKM','3A4P', '3CCN', '3CD8','3CE3', '3CTH', '3CTJ', '3DKC', '3DKF', '3DKG', '3EFJ', '3EFK', '3F66', '3F82', '3I5N','3L8V', '3LQ8', '3Q6U', '3Q6W', '3QTI', '3R7O', '3RHK', '3U6H', '3U6I', '3VW8', '3ZBX', '3ZC5', '3ZCL', '3ZXZ', '3ZZE', '4AOI', '4AP7', '4DEG', '4DEH', '4DEI', '4EEV', '4GG5', '4GG7', '4IWD', '4KNB', '4MXC', '4R1V', '4R1Y', '4XMO', '4XYF', '5DG5', '5EOB', '5EYC', '5EYD', '5HLW', '5HNI', '5HO6', '5HOA', '5HOR', '5HTI', '5T3Q', '5UAB', '5UAD', '5UAF', '5YA5', '6SD9', '6SDC', '6SDD', '6SDE', '6UBW', '7B3Z', '7B3Q', '7B42', '7B41', '7B44', '7B43', '7B3T', '7B3W', '7B3V', '7B40', '7V3R', '7V3S', '7Y4T', '7Y4U', '8AN8','8ANS', '8OUU', '8OUV', '8OU7', '8OVZ', '8OW3', '8OWG']

L = []
i=0

while i<100: #len(pdbs):
  #create a conditional statement in case structure is missing --> avoid error for 'NoneType' object
  b = resi_1066_xyz_dictionary.get(pdbs[i])                                     # assigns specific coordinate array for residue to variable
  c = resi_1129_xyz_dictionary.get(pdbs[i])
  d = resi_1062_xyz_dictionary.get(pdbs[i])
  e = resi_1125_xyz_dictionary.get(pdbs[i])
  f = resi_1058_xyz_dictionary.get(pdbs[i])
  g = resi_1121_xyz_dictionary.get(pdbs[i])

  #l = resi_1124_xyz_dictionary.get(pdbs[i])
  #m = resi_1153_xyz_dictionary.get(pdbs[i])
  #n = resi_1112_xyz_dictionary.get(pdbs[i])
  #print (c)
  #print (b)
  #print (c[0])
  #print (b[0])

  if b and c and d and e and f and g is not None:
    dx = d[0]-e[0]
    dy = d[1]-e[1]
    dz = d[2]-e[2]
    h = math.sqrt(dx**2 + dy**2 + dz**2)
    print (h)
    ddx = b[0]-c[0]
    ddy = b[1]-c[1]
    ddz = b[2]-c[2]
    hh = math.sqrt(ddx**2 + ddy**2 + ddz**2)
    print (hh)
    dddx = f[0]-g[0]
    dddy = f[1]-g[1]
    dddz = f[2]-g[2]
    hhh = math.sqrt(dddx**2 + dddy**2 + dddz**2)
    print (hhh)
    #line = re.findall(r'resolution', pdbs[i]
    dis = [pdbs[i],h,hh,hhh]
    L.append(dis)
    print (i)
    i = i+1
  else:
    i = i+1

with open('JM_MET_aligned_strutures_2.csv', 'w', newline='') as file:
  writer = csv.writer(file)
  writer.writerow(['PDB','1062 to 1125 distance', '1066 to 1129 distance', '1058 to 1121 distance'])
  writer.writerows(L)

python end


print(L)


#distance i. 1066 and n. CA, i. 1129 and n. CA
#distance i. 1062 and n. CA, i. 1125 and n. CA
#distance i. 1058 and n. CA, i. 1121 and n. CA


set_view (\
     0.905076325,    0.277848482,    0.321927637,\
    -0.020534338,   -0.727594256,    0.685700357,\
     0.424753398,   -0.627221465,   -0.652822793,\
    -0.000003025,    0.000030279, -299.019836426,\
    13.646133423,   10.031967163,  151.169708252,\
    79.283798218,  518.755126953,  -20.000000000 )
