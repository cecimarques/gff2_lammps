# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 13:44:53 2025

@author: u2371853
"""

import numpy as np

# INPUT: CHANGE ACCORDINGLY 
# Please put here the number of molecules in the simulation domain and 
# the number of atoms in a molecule
number_of_molecules = 343
NA = 14

# INPUT: CHANGE ACCORDINGLY 
# The file dihedralcoeff is a file you should build manually by cp'ing the
# dihedral potential parameters output in the main code. Below you should
# put the number of lines in this file (which should be >= the number of di-
# hedrals in a molecule)
number_of_dihedrals_coeff = 33

dihedrals_coeff = []
ofi = open("dihedralcoeff", 'r')
for it_1 in range(0, number_of_dihedrals_coeff):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(4)):
            dihedrals_coeff.append(float(e))
dihedrals_coeff = np.array(dihedrals_coeff,float)
dihedrals_coeff = dihedrals_coeff.reshape(number_of_dihedrals_coeff,4)

# INPUT: CHANGE ACCORDINGLY 
# Please put here the number of dihedrals in the 1-molecule LigParGen file.
# This section, with the column separted by tab bars, compose the file "di-
# hedrals" below.
number_of_dihedrals = 27

dihedrals = []
ofi = open("dihedrals", 'r')
for it_1 in range(0, number_of_dihedrals):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(6)):
            dihedrals.append(float(e))
dihedrals = np.array(dihedrals,float)
dihedrals = dihedrals.reshape(number_of_dihedrals,6)

# --------------------------------------------------------------------
# First, lets see which of these dihedrals are printing twice the same potential
# for a given dihedral type simply as a consequence of reading the same line of
# the ffdihedral file twice (for which in the case of bonds and angles I deleted
# repetitions manually).

# The array line_to_delete will take the lines of the dihedralcoeff file that 
# I need to delete as they are formed by dihedrals of the given type that has 
# the exact same parameters appearing in the file, case which I am going to as-
# sume to be caused by the dihedral being a palindrome (afterall, why set two 
# potentials with the exact same parameters instead of simply multiplying by 2 
# the force constant of the other one?). 
line_to_delete = np.array([[0]])
for it_1 in range(0, len(dihedrals_coeff)):
    given_type = dihedrals_coeff[it_1,0]
    for it_2 in range(it_1 + 1, len(dihedrals_coeff)):
        if (dihedrals_coeff[it_2,0] == dihedrals_coeff[it_1,0]) & (dihedrals_coeff[it_2,1] == dihedrals_coeff[it_1,1]) & (dihedrals_coeff[it_2,2] == dihedrals_coeff[it_1,2]) & (dihedrals_coeff[it_2,3] == dihedrals_coeff[it_1,3]):
            # ---------------------------------------------------------
            # Before accumulating the line I am supposed to delete, I will check
            # if I havent already stored it.
            flag = 0
            for it_3 in range(0, len(line_to_delete)):
                if line_to_delete[it_3,0] == it_2:
                    flag = 1
                    # No need to continue the it_3 for loop in this case:
                    break
            # ---------------------------------------------------------
            if flag == 0:
               line_to_delete = np.concatenate((line_to_delete, np.array([[it_2]])), axis = 0)

# Lets delete the line containing the 0 element.
line_to_delete = np.delete(line_to_delete, 0, 0)

# Now lets delete the lines marked in the line_to_delete array.
# My methodology for deletion requires first ordering the line_to_deliete
# array in descending order.
for it_1 in range(0, len(line_to_delete)):
    for it_2 in range(it_1 + 1, len(line_to_delete)):
        if line_to_delete[it_1,0] > line_to_delete[it_2,0]:
            temporary = line_to_delete[it_1,0]
            line_to_delete[it_1,0] = line_to_delete[it_2,0]
            line_to_delete[it_2,0] = temporary
for it_1 in reversed(range(0, len(line_to_delete))):
    given_line = int(line_to_delete[it_1,0])
    dihedrals_coeff = np.delete(dihedrals_coeff, given_line, 0)
    
# --------------------------------------------------------------------------------
# Now lets read again the dihedral_coeff array and store information concerning addi-
# tional dihedral potentials, and thus new types, that need to be tuned.
# The array duplicate will keep the dihedral type index that has multiple rightful 
# declarations in the dihedral_coeff array and the line containing the parameters of
# this additional potential to-be-tuned.
# Ultimately, each line in the duplicate array implies an additional dihedral poten-
# tial that needs to be tuned for the dihedral of the given corresponding type.
duplicate = np.zeros((1,2))
for it_1 in range(0, len(dihedrals_coeff)):
    given_type = dihedrals_coeff[it_1,0]
    for it_2 in range(it_1 + 1, len(dihedrals_coeff)):
        if dihedrals_coeff[it_2,0] == given_type:
            if (dihedrals_coeff[it_1,1] != dihedrals_coeff[it_2,1]) or (dihedrals_coeff[it_1,2] != dihedrals_coeff[it_2,2]) or (dihedrals_coeff[it_1,3] != dihedrals_coeff[it_2,3]):
                # If I found a line having another potential for this given dihedral type,
                # I will store the index of this line (and the type of the dihedral). 
                # --------------------------------------------------------
                # But note that I may have already identified this line before, so I 
                # will check if indeed this is the case before storing info on it a-
                # gain.
                flag = 0
                for it_3 in range(0, len(duplicate)):
                    if duplicate[it_3,1] == it_2:
                        flag = 1
                # --------------------------------------------------------
                if flag == 0:
                    duplicate = np.concatenate((duplicate, np.array([[given_type, it_2]])), axis = 0)

# Deleting the initial line used to initialize the array:
duplicate = np.delete(duplicate, 0, 0)

# Now I am going to create new dihedral types for each additional potential that I 
# need to tune over a given dihedral.
# First I am concatenating an "empty" column in the duplicate array that will store
# the ID of the new type.
duplicate = np.concatenate((duplicate, np.zeros((len(duplicate),1))), axis = 1)
dihedral_type = len(dihedrals) + 1
for it_1 in range(0, len(duplicate)):
    given_line = int(duplicate[it_1,1])
    dihedrals_coeff[given_line,0] = dihedral_type
    duplicate[it_1,2] = dihedral_type
    dihedral_type = dihedral_type + 1

# Now I need to duplicate the dihedrals in the dihedral section and assign them the
# new type. I am doing this for the 1-molecule dihedral section as obtained from
# LigParGen. 
dihedral_ID = len(dihedrals) + 1
for it_1 in range(0, len(duplicate)):
    given_type = duplicate[it_1,0]
    for it_2 in range(0, len(dihedrals)):
        if dihedrals[it_2,0] == given_type:
            temporary_array = np.zeros((1,6))
            temporary_array[0,:] = dihedrals[it_2,:]
            temporary_array[0,0] = dihedral_ID
            temporary_array[0,1] = duplicate[it_1,2]
            dihedrals = np.concatenate((dihedrals, temporary_array), axis = 0)
            dihedral_ID = dihedral_ID + 1
            break

# Now lets "propagate" the dihedral section for the given number of molecules I
# have in the simulation domain. The setup I have below for declaring the atom 
# IDs in the "propagated" dihedral section works well because of the way I have
# designed my code to propagate the 1-molecule template across space.
number_of_dihedrals = len(dihedrals)
non_propagated_dihedral_length = len(dihedrals)
for it_1 in range (1, number_of_molecules):
    for it_2 in range (0, number_of_dihedrals):
        temporary_array = np.zeros((1,6))
        temporary_array[0,0] = (number_of_dihedrals)*it_1 + dihedrals[it_2,0]
        temporary_array[0,1] = dihedrals[it_2,1]
        temporary_array[0,2] = NA*it_1 + dihedrals[it_2,2]
        temporary_array[0,3] = NA*it_1 + dihedrals[it_2,3]
        temporary_array[0,4] = NA*it_1 + dihedrals[it_2,4]
        temporary_array[0,5] = NA*it_1 + dihedrals[it_2,5]
        dihedrals = np.concatenate((dihedrals, temporary_array), axis = 0)

# -------------------------------------------------------------------------------
# This bit below counts the number of dihedrals types that do not appear in the Dihedral
# Coeff section, which should be a consequence of parameters for dihedral potentials for-
# med by the given underlying sequence of CGenFF atom types not existing in the ffdihedral
# file *OR* being declared with a X for some of the atom types forming it to indicate the
# potential parameters do not depend on them. In the latter, I will later need to check
# the ffdihedral file later and in the former, I will need to "tune" a potential there and
# put a 0 in the force constant.

# The variable counter will count the number of missing dihedrals
counter = 0
# This list will take the index of dihedral types that do not appear declared in the 
# Dihedral Coeff section as I find them.
missing_type = []

for it_1 in range(0, non_propagated_dihedral_length):
    given_type = it_1 + 1 
    flag = 0
    for it_2 in range(0, len(dihedrals_coeff)):
        if dihedrals_coeff[it_2,0] == given_type:
            flag = 1
    # if I never found it a dihedral of the given type, I need to take action:
    if flag == 0:
        counter = counter + 1
        missing_type.append(float(given_type))

print(missing_type)

# This will add a line in the Dihedral Coeff section with the given dihedral
# type ID and parameters that lead to effectively no dihedral being tuned.
# Note that you need to change that manually later in case you end up reali-
# zing there is a dihedral potential in the ffdihedral file with a "X" decla-
# red instead of an atom type.
temporary_array = np.zeros((counter,4))
for it_1 in range(0, counter):
    temporary_array[it_1,0] = missing_type[it_1]
dihedrals_coeff = np.concatenate((dihedrals_coeff, temporary_array), axis = 0)

# ----------------------------------------------------------------------------
# Now I need to tune in the value of the w parameter in the dihedral_style used in
# LAMMPS, which will be used to tune in the 1-4 interactions. These should be 1 for
# all dihedral potentials except for those that act on the same dihedrals OR those 
# that act on a different set of atoms but the 1 and 4 atoms forming the dihedral 
# are the same as another dihedral potential.

# Beforehand, I will need to order the dihedrals_coeff in ascending order of dihe-
# dral type in order for the methodology below to work well.
temporary_array = np.zeros((1,4))
for it_1 in range(0, len(dihedrals_coeff)):
    for it_2 in range(it_1 + 1, len(dihedrals_coeff)):
        if dihedrals_coeff[it_1,0] > dihedrals_coeff[it_2,0]:
            temporary_array[0,:] = dihedrals_coeff[it_1,:]
            dihedrals_coeff[it_1,:] = dihedrals_coeff[it_2,:]
            dihedrals_coeff[it_2,:] = temporary_array[0,:]
            
# This column will hold the value of w for each dihedral potential, and initially I
# will set all of them to 1. I will subsequentially 0 those that are needed.
dihedrals_coeff = np.concatenate((dihedrals_coeff, np.zeros((len(dihedrals_coeff),1))), axis = 1)
dihedrals_coeff[:,4] = 1
for it_1 in range(0, non_propagated_dihedral_length):
    atom1 = dihedrals[it_1,2]
    atom4 = dihedrals[it_1,5]
    for it_2 in range(it_1 + 1, non_propagated_dihedral_length):
        # -------------------------------------------------------------
        if (dihedrals[it_2,2] == atom1) & (dihedrals[it_2,5] == atom4):
            type_index = dihedrals[it_2,1]
            dihedrals_coeff[int(type_index - 1),4] = 0
        # -------------------------------------------------------------
        if (dihedrals[it_2,2] == atom4) & (dihedrals[it_2,5] == atom1):
            type_index = dihedrals[it_2,1]
            dihedrals_coeff[int(type_index - 1),4] = 0


