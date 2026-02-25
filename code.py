# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:46:57 2025

@author: u2371853
"""

import numpy as np

atomtype = []
# INPUT: CHANGE ACCORDINGLY 
# Here you must put the CGenFF atom types for all atoms composing your mo-
# lecule: you can find this in the gmx.top file you get from 
# https://cgenff.com/.
# These atom types should be informed in ascending order of the LAMMPS a-
# tom ID in the 1-molecule LAMMPS data file got from LigParGen. 
for it_1 in range(6):
    atomtype.append("CG2R61")
for it_1 in range(6):
    atomtype.append("HGR61")

NA = len(atomtype)

# INPUT: CHANGE ACCORDINGLY
# Here you must put the charges for each of the atom types composing your
# molecule. Just like before, you should find this in the gmx.top file you
# get from https://cgenff.com/.
# The charges should be informed in ascending order of the LAMMPS atom ID 
# in the 1-molecule LAMMPS data file got from LigParGen. 
charge = [-0.115, -0.115, -0.115, -0.115, -0.115, -0.115, 0.115, 0.115, 0.115, 0.115, 0.115, 0.115]

# INPUT : CHANGE ACCORDINGLY
# Here are further things you need to change. The number of bonds, angles, 
# dihedrals and improper found declared in the 1-molecule LAMMPS data file
# you get from LigParGen. The number_of_molecules is the number of molecu-
# les you have in *another* LAMMPS data file, which you will actually use
# in the simulation, that contains a given number of the 1-molecule pro-
# pagated in space.
number_of_molecules = int(8*8*8)
number_of_bonds = 12
number_of_angles = 18
number_of_dihedrals = 24
number_of_impropers = 6
# ---------------------------------------------------------------------------------
epsilon = []
epsilon14 = []
sigma = []
sigma14 = []

# This is the number of lines in the ffnonbonded file (see github READ.me 
# for the meaning of this file).
# I believe there is no need to change this, as the file cgenff output for
# no matter what compound should be the same.
lines_ffnonbonded = 458

# Now I am going to find which lines of the ffnonbonded.itp file I care about wi-
# thin the scope of this given compound AND save the "atom name", epsilon and 
# sigma that I find there for each of the atom IDs. Later I will be using this
# to generate the LJ parameters via mixing rules.
for it_1 in range(0, len(atomtype)):
    value = atomtype[it_1]
    # ----------------------------------------------------------------------------
    ofi = open("ffnonbonded", 'r')
    for it_2 in range(0, lines_ffnonbonded):
        dump = ofi.readline()
        # This variable will count how many non-zero string sequences (i.e. text) in 
        # the dump variable I will find.
        count = 0
        # This variable will take value of 0 everytime I read a space bar or \t in the dump
        # file as well as will be 0 for the first time I read a non-space-bar string slot: 
        # in that latter case, the value of flag will be 1 for all non-space-bar slot that
        # I read in the sequence.
        flag = 0
        # This variable functions to me as an indicator on whether or not I should be saving
        # data belonging to a text in the dump. Naturally it should only be "activated" if I in-
        # deed met a line of the ffnonbonded.itp file that has the atom type corresponding to
        # the one that appears in the line it_1 of the list atomtype
        save_data = 0
        # This loops over all slots of the dump string
        for it_3 in range(0,len(dump)):
            if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                flag = 0
            if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                # ------------------------------------------------------------
                if (flag == 0):
                    # Found a text in the line!
                    count = count + 1
                    # defining an empty "string" variable.
                    value = str()
                    flag = 1
                # ------------------------------------------------------------
                # This concatenates the element in the slot it_3 of dump to the variable
                # "value".
                value = value+''.join(dump[it_3])
                
                # This runs everytime I am *not* in the last iteration of it_3.
                if it_3 != (len(dump)-1):
                    # If the next slot of the dump is a space bar or a tab, for sure I will be in the 
                    # last slot that belongs to this text.
                    # I will check if it is a text that matches what is the atomtype and possibly set 
                    # the save_data value accordingly. The value of save_data will be 1 once I get to
                    # the values of epsilon and sigma in the file.
                    # The if condition above was necessary because in the last slot of dump I will 
                    # have a problem in the if condition below (out of range): I will be taking ca-
                    #re of that scenario below.
                    if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t') or (dump[it_3+1] == '\n'):
                        if value == atomtype[it_1]:
                            save_data = 1
                        if (save_data == 1) & (count == 6):
                            sigma.append(value)
                        if (save_data == 1) & (count == 7):
                            epsilon.append(value)
                            # lets reset the save_data value since I already saved all
                            # the data I needed to save and break the loop of it_3 as
                            # I no longer will need to run it for the atom correspondnig
                            # to atomtype[it_1]
                            save_data = 0
                            break
                        
                        
                # The only possible scenario of having it_3 == (len(dump)-1) in which I need to do 
                # something, is if the remainings of the decimals of epsilon are in the last slot. 
                # The condition I have set should work as save_data will only still be 1 if the la-
                # st character in the line is *not* a space bar or a tab or an enter BUT actually
                # is a random character, such as ; for example, that is "glued" to the 8th text (i.
                # e. epsilon value) that appears in the line. This should not happen, as an enter
                # will always occur, and so this line should never yield an append.
                if it_3 == (len(dump)-1):
                    if (save_data == 1) & (count == 7):
                        epsilon.append(value)
                        save_data = 0
                        break

# ----------------------------------------------------------------------
# Now lets find the parameters for the 1-4 interaction frmo the [ pairtypes ] directi-
# ve of the GROMACS input file which is cp'd into a new file named ffpairtype. 
# I am hoping that the parameters for 1-4 interaction between all atom types is listed
# here, and if there is something missing (e.g. if generated via the gen-yes keyword in
# gromacs), I will also account for the need to generate it via the mixing rules in case 
# I dont find it for a given pair.
# Furthermore, I want to tell you that I am going to account for repetition of para-
# meters for palindrome pair of atom types (which should, naturally, happen when the
# two atom types forming the pair is the same).

lines_ffnpairtype = 41769

for it_1 in range(0, len(atomtype)):
    ID1 = atomtype[it_1]
    for it_1i in range(it_1, len(atomtype)):
        ID2 = atomtype[it_1i]
        # --------------------------------------------------------------
        # This variable will mark to me when I find and when I dont find the parameters
        # for the given pair of atom types in the list of [ pairtypes ]. If I dont find
        # it, I will generate it.
        ever_ran = 0
        # --------------------------------------------------------------
        ofi = open("ffpairtype", 'r')
        for it_2 in range(0, lines_ffnpairtype):
            dump = ofi.readline()
            count = 0
            flag = 0
            # These variables will later take the value of the atom names that appear
            # in a given line of the ffpairtype file.
            ID1_fileff = 0
            ID2_fileff = 0
            for it_3 in range(0,len(dump)):
                if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                    flag = 0
                if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                    # ------------------------------------------------------------
                    if (flag == 0):
                        count = count + 1
                        value = str()
                        flag = 1
                    # ------------------------------------------------------------
                    value = value+''.join(dump[it_3])
                    
                    if it_3 != (len(dump)-1):
                        if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t'):
                            if count == 1:
                                ID1_fileff = value
                            if count == 2:
                                ID2_fileff = value
                                # No need to continue reading this line if I finally found
                                # the second ID declared in the bond in the ffbond file.
                                break
            # This condition will only occur if I happened to have broken the loop in the
            # context of the previous command line (no other possibility), which should o-
            # curr for all lines that I read in the ffpairtype file.
            if count == 2:
                # The if conditions below that will only enter if I happen to have read the line 
                # of the ffpairtype file that contains the information for the two given "atom na-
                # mes". In that latter case I will make sure I output the information already in 
                # the LAMMPS format.
                if (ID1 == ID1_fileff) & (ID2 == ID2_fileff) & (ever_ran == 0):
                    # ------------------------------------------------------
                    # This pair has coefficients listed.
                    ever_ran = 1
                    # ------------------------------------------------------
                    # These for now empty strings will later take the value of
                    # the parameters.
                    value1 = str()
                    value2 = str()
                    # This variable will count the text in the line.
                    count2 = 0
                    flag = 0
                    for it_4 in range(0, len(dump)):
                        if (dump[it_4] == ' ') or (dump[it_4] == '\t') or (dump[it_4] == '\n') or (dump[it_4] == '') or (dump[it_4] == ""):
                            flag = 0
                            continue
                        # I should only get to the command line below if I do not
                        # fall into the if condition above.
                        # --------------------------------------------------
                        # If the condition below is met, I am in the first charac-
                        # ter of the text.
                        if flag == 0:
                            count2 = count2 + 1
                            flag = 1
                        # ---------------------------------------------------
                        # These lines will store all the characters of the text
                        # that appear in the 5th and 4th columns of the file.
                        if (count2 == 5) & (flag == 1):
                            value1 = value1+dump[it_4]
                        if (count2 == 4) & (flag == 1):
                            value2 = value2+dump[it_4]
                    # ---------------------------------------
                    epsilon14.append(value1)
                    sigma14.append(value2)
                    # ---------------------------------------
                # This accounts for the other possibility of atomtypes.
                if (ID2 == ID1_fileff) & (ID1 == ID2_fileff) & (ever_ran == 0):
                    # ------------------------------------------------------
                    ever_ran = 1
                    # ------------------------------------------------------
                    value1 = str()
                    value2 = str()
                    count2 = 0
                    flag = 0
                    for it_4 in range(0, len(dump)):
                        if (dump[it_4] == ' ') or (dump[it_4] == '\t') or (dump[it_4] == '\n') or (dump[it_4] == '') or (dump[it_4] == ""):
                            flag = 0
                            continue
                        # --------------------------------------------------
                        if flag == 0:
                            count2 = count2 + 1
                            flag = 1
                        # ---------------------------------------------------
                        if (count2 == 5) & (flag == 1):
                            value1 = value1+dump[it_4]
                        if (count2 == 4) & (flag == 1):
                            value2 = value2+dump[it_4]
                    # ---------------------------------------
                    epsilon14.append(value1)
                    sigma14.append(value2)
                    # ---------------------------------------
        # Lets generate the coefficients for this 1-4 pair if I finished
        # reading the ffpairtype file AND did not find a match for them
        # in any line.
        if ever_ran == 0:
            value1 = (float(epsilon[it_1])*float(epsilon[it_1i]))**(1/2)
            value2 = (float(sigma[it_1]) + float(sigma[it_1i]))/2
            epsilon14.append(value1)
            sigma14.append(value2)
                    
# -----------------------------------------------------------------------------
# Now lets print the converted parameter values for each atom type that exists in
# the file.

epsilon = np.array(epsilon,float)
epsilon = epsilon.reshape(NA,1)
sigma = np.array(sigma,float)
sigma = sigma.reshape(NA,1)

# Unit and GROMACS->LAMMPS potential form conversion:
epsilon = epsilon/4.184
sigma = sigma*10

# This will count how many lines the PairIJ section should have in LAMMPS
shape = 0
for it_1 in reversed(range(1, len(atomtype) + 1)):
    shape = shape + it_1
shape = int(shape)
output_epsilon = np.zeros((shape,1))
output_sigma = np.zeros((shape,1))

# Here I am generating the parameters for pairwise interactions other than 
# 1-4 via mixing rules.
given_line = 0
for it_1 in range(0, len(atomtype)):
    for it_2 in range(it_1, len(atomtype)):
        output_epsilon[given_line,0] = (epsilon[it_1,0]*epsilon[it_2,0])**(1/2)
        output_sigma[given_line,0] = (sigma[it_1,0] + sigma[it_2,0])/2
        given_line = given_line + 1

# Now I need to do the processing of the epsilon and sigma that are specific for the
# 1-4 interactions. I am expecting to have already have *all* the coefficients when
# scanning the pairlist: if the amount of coefficients is less or, for some reason,
# more than expected, there will be an error printed below as it will hopefully be 
# impossible to generate an array with the dimensions I am setting.
epsilon14 = np.array(epsilon14,float)
epsilon14 = epsilon14.reshape(shape,1)
sigma14 = np.array(sigma14,float)
sigma14 = sigma14.reshape(shape,1)

# Unit and GROMACS->LAMMPS potential form conversion:
epsilon14 = epsilon14/4.184
sigma14 = sigma14*10

# Now I am printing these non-bonded parameters for me to directly copy and paste to
# the LAMMPS data file.
nb_output = np.concatenate((output_epsilon, output_sigma, epsilon14, sigma14), axis = 1)
prefix = np.zeros((len(nb_output),2))
line = 0
for it_1 in range(0, len(atomtype)):
    for it_2 in range(it_1, len(atomtype)):
        prefix[line,0] = it_1 + 1
        prefix[line,1] = it_2 + 1
        line = line + 1
print("non-bonded parameters")
for it_1 in range(0, len(nb_output)):
    print(int(prefix[it_1,0]),"\t", int(prefix[it_1,1]),"\t", "{0:.6f}".format(nb_output[it_1,0]),"\t", "{0:.6f}".format(nb_output[it_1,1]), "\t", "{0:.6f}".format(nb_output[it_1,2]), "\t", "{0:.6f}".format(nb_output[it_1,3]))

# -------------------------------------------------------------------------------
# Now the processing of the charges: the goal here is just to write a file
# containing in 1 column the charges*number_of_molecules to copy paste into
# the LAMMPS data file.

# Lets transform it into an array:
charge = np.array(charge,float)
charge = charge.reshape(NA,1)

# Sanity check to see if my charges are summing to one. Please check the output in the
# kernel.
charsum = 0
for it_1 in range(0, len(charge)):
    charsum = charsum + charge[it_1]
print("Molecule net charge:", charsum)

# Now let's output the 1-column file containing the charges. Note that the
# approach for declaring the charge of each atom set up below works well be-
# cause I propagated the molecule in space by having the IDs of a given atom
# be spaced by NA (i.e. number of atom in a molecule) from its ghost.
output_charge = np.zeros((1,1))
for it_1 in range(0, number_of_molecules):
    output_charge = np.concatenate((output_charge, charge), axis = 0)
output_charge = np.delete(output_charge, 0, 0)
ofi = open("chargesff.out", 'w')   
for it_1 in range(len(output_charge)):
    for it_2 in range(0,1):
            ofi.write(str("{0:.6f}".format(output_charge[it_1,it_2])))
    ofi.write('\n')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# This reads the bond section of the LigParGen file except that the spacing bet-
# ween the columns is design to be a \t.
bonds = []
ofi = open("bonds", 'r')
for it_1 in range(0, number_of_bonds):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(4)):
            bonds.append(float(e))
bonds = np.array(bonds,float)
bonds = bonds.reshape(number_of_bonds,4)

# INPUT: CHANGE ACCORDINGLY
# This is the number of lines under the [ bondtypes ] directive of the file 
# ffbonded.itp file.
# You should change this accordinglyas you should take into consideration not
# not only the bond types in the standard ffbonded.itp file (which should be
# the same for all compounds) but also you need to add the bond types that mi- 
# ght be found in the file dra_ffbonded.itp for ur specific molecule.
lines_ffbond = 1451

print("\n")
print("Bond Potential Parameters")
for it_1 in range(0, len(bonds)):
    # These receive the atom types of the two atoms declared as bonded in 
    # the line it_1 of the file "bonds".
    ID1 = atomtype[int(bonds[it_1,2]-1)]
    ID2 = atomtype[int(bonds[it_1,3]-1)]
    # ----------------------------------------------------------------------------
    ofi = open("ffbond", 'r')
    for it_2 in range(0, lines_ffbond):
        dump = ofi.readline()
        count = 0
        flag = 0
        # These variables will later take the value of the atom types that appear
        # in a given line of the ffbond file.
        ID1_fileff = 0
        ID2_fileff = 0
        for it_3 in range(0,len(dump)):
            if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                flag = 0
            if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                # ------------------------------------------------------------
                if (flag == 0):
                    count = count + 1
                    value = str()
                    flag = 1
                # ------------------------------------------------------------
                value = value+''.join(dump[it_3])
                
                if it_3 != (len(dump)-1):
                    if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t'):
                        if count == 1:
                            ID1_fileff = value
                        if count == 2:
                            ID2_fileff = value
                            # No need to continue reading this line of the ffbonds
                            # file if I finally found the second ID declared in the 
                            # bond in the ffbond file.
                            break
        # This condition will only occur if I happened to have broken the loop in the
        # context of the previous command line (no other possibility), which should o-
        # curr for all lines that I read in the ffbond file though. 
        if count == 2:
            # The if conditions below that will only enter if I happen to have read the line 
            # of the ffbond file that contains the information on the bond between the two 
            # given "atom names"). In that latter case I will make sure I output the informa-
            # tion already in the LAMMPS format.
            if (ID1 == ID1_fileff) & (ID2 == ID2_fileff):
                # These for now empty strings will later take the value of
                # the bond potential parameters.
                bit1 = str()
                bit2 = str()
                # This variable will count the text in the line.
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t') or (dump[it_4] == '\n'):
                        flag = 0
                        continue
                    # I should only get to the command linea below if I do not
                    # fall into the if condition above.
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    # These lines will store all the characters of the text
                    # that appear in the 5th and 4th text slot of the line.
                    if (count2 == 5) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 4) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                # ---------------------------------------
                # Converting the units and accounting for differences in poten-
                # tial form.
                bit1 = float(bit1)
                bit1 = "{0:.6f}".format(bit1/(4.184*100*2))
                bit1 = str(bit1)
                bit2 = float(bit2)
                bit2 = "{0:.6f}".format(bit2*10)
                bit2 = str(bit2)
                # ---------------------------------------
                out = str(int(bonds[it_1,1]))+str("\t")+bit1+str("\t")+bit2
                print(out)
            # Now the same as above for the possibility of inverted IDs.
            if (ID2 == ID1_fileff) & (ID1 == ID2_fileff):
                bit1 = str()
                bit2 = str()
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t') or (dump[it_4] == '\n'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 5) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 4) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                # ---------------------------------------
                bit1 = float(bit1)
                bit1 = "{0:.6f}".format(bit1/(4.184*100*2))
                bit1 = str(bit1)
                bit2 = float(bit2)
                bit2 = "{0:.6f}".format(bit2*10)
                bit2 = str(bit2)
                # ---------------------------------------
                out = str(int(bonds[it_1,1]))+str("\t")+bit1+str("\t")+bit2
                print(out)

# There will be repeated occurances if these occur in the .itp files OR if
# the atom names of the two IDs are the same. You will be able to spot re-
# petition by seeing that more than one set of parameter is printed for a sa-
# me bond type. If the parameters are the same, I suppose it is no problem
# since it implies that the info through the file for the given "atom names" 
# is coherent. If not, it means there is a problem and you will need to go 
# check the ffbond file to see if there is some comment on the lines that 
# help you choose which set of parameters if more suitable for your system.
# Note also that obviously you will need to erase the repetition after copy-
# ing pasting everything to the LAMMPS file. This is great, because it will
# then force you to check.

# ------------------------------------------------------------------------------
# The idea for this bit is the same as for the bonds above, so I wont repeat
# the comments.
angles = []
ofi = open("angles", 'r')
for it_1 in range(0, number_of_angles):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(5)):
            angles.append(float(e))
angles = np.array(angles,float)
angles = angles.reshape(number_of_angles,5)

lines_ffangle = 5043

print("\n")
print("Angle Potential Parameters")
for it_1 in range(0, len(angles)):
    ID1 = atomtype[int(angles[it_1,2]-1)]
    ID2 = atomtype[int(angles[it_1,3]-1)]
    ID3 = atomtype[int(angles[it_1,4]-1)]
    # ----------------------------------------------------------------------------
    ofi = open("ffangle", 'r')
    for it_2 in range(0, lines_ffangle):
        dump = ofi.readline()
        count = 0
        flag = 0
        ID1_fileff = 0
        ID2_fileff = 0
        ID3_fileff = 0
        for it_3 in range(0,len(dump)):
            if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                flag = 0
            if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                # ------------------------------------------------------------
                if (flag == 0):
                    count = count + 1
                    value = str()
                    flag = 1
                # ------------------------------------------------------------
                value = value+''.join(dump[it_3])
                
                if it_3 != (len(dump)-1):
                    if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t'):
                        if count == 1:
                            ID1_fileff = value
                        if count == 2:
                            ID2_fileff = value
                        if count == 3:
                            ID3_fileff = value
                            break

        if count == 3:
            if (ID1 == ID1_fileff) & (ID2 == ID2_fileff) & (ID3 == ID3_fileff):
                bit1 = str()
                bit2 = str()
                bit3 = str()
                bit4 = str()
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 6) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 5) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                    if (count2 == 8) & (flag == 1):
                        bit3 = bit3+dump[it_4]
                    if (count2 == 7) & (flag == 1):
                        bit4 = bit4+dump[it_4]
                # ---------------------------------------
                bit1 = float(bit1)
                bit1 = "{0:.6f}".format(bit1/(4.184*2))
                bit1 = str(bit1)
                bit2 = float(bit2)
                bit2 = "{0:.6f}".format(bit2)
                bit2 = str(bit2)
                bit3 = float(bit3)
                bit3 = "{0:.6f}".format(bit3/(4.184*100*2))
                bit3 = str(bit3)
                bit4 = float(bit4)
                bit4 = "{0:.6f}".format(bit4*10)
                bit4 = str(bit4)
                # ---------------------------------------
                out = str(int(angles[it_1,1]))+str("\t")+bit1+str("\t")+bit2+str("\t")+bit3+str("\t")+bit4
                print(out)
            
            # Taking into account the swap of atom names in the 1st and 3rd
            # text slot:
            if (ID3 == ID1_fileff) & (ID2 == ID2_fileff) & (ID1 == ID3_fileff):
                bit1 = str()
                bit2 = str()
                bit3 = str()
                bit4 = str()
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 6) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 5) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                    if (count2 == 8) & (flag == 1):
                        bit3 = bit3+dump[it_4]
                    if (count2 == 7) & (flag == 1):
                        bit4 = bit4+dump[it_4]
                # ---------------------------------------
                bit1 = float(bit1)
                bit1 = "{0:.6f}".format(bit1/(4.184*2))
                bit1 = str(bit1)
                bit2 = float(bit2)
                bit2 = "{0:.6f}".format(bit2)
                bit2 = str(bit2)
                bit3 = float(bit3)
                bit3 = "{0:.6f}".format(bit3/(4.184*100*2))
                bit3 = str(bit3)
                bit4 = float(bit4)
                bit4 = "{0:.6f}".format(bit4*10)
                bit4 = str(bit4)
                # ---------------------------------------
                out = str(int(angles[it_1,1]))+str("\t")+bit1+str("\t")+bit2+str("\t")+bit3+str("\t")+bit4
                print(out)

# Same comment found in the end of the section for finding the bond potential para-
# meters hold here for the angles.

# ------------------------------------------------------------------------------
# The idea for this bit is the same as for the bonds above, so I wont repeat
# the comments. There are however two differences here: the first and main
# one is that CGenFF prescribes sometimes more than a single dihedral poten-
# tial of the form given in https://docs.lammps.org/dihedral_charmm.html for
# a single dihedral. In other words, this is not some problem related to pa-
# lindrome sequences or anything. So you should basically copy paste the out-
# put you get for dihedral potentials, without making any modification, in-
# to a file called dihedralcoeff, which will serve as input for a second co-
# de named dihedral_solution_CGenFF.py (see maintree page on github). This
# leads to the second difference previously mentioned: since there is a co-
# de devoted to tackling this particularity of CGenFF, I also made sure the
# code takes care of deleting palindrome-repeated dihedral potentials, whi-
# ch are believed to be a result of the problem mentioned above for the
# bonds instead of 2 different dihedral potentials being tuned (otherwise
# it would be necessary to simply use 2x the value of K instead of tunning
# two different dihedral potentails of the dihedral_style charmm).

# Also I would like to note that I have not build this part of the code to
# address dihedral potentails tuned over dihedrals where one or more atom
# type in the file is declared as "X" as a consequence of it not being im-
# portant for the parametrization. As this is the case for some dihedrals
# within CGenFF (see the ffdihedral file), there will be missing dihedral
# types: this will also be evidenced in the separate code I wrote to tack-
# le the preivously mentioned problem.
# Finally, note that the last dihedral potential parameter for the dihed-
# dral_style charmm will still be missing as it shouldnt exist in the .itp
# file: this will be taken care of also in the dihedral_solution_CGenFF.py
# code.

dihedrals = []
ofi = open("dihedrals", 'r')
for it_1 in range(0, number_of_dihedrals):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(6)):
            dihedrals.append(float(e))
dihedrals = np.array(dihedrals,float)
dihedrals = dihedrals.reshape(number_of_dihedrals,6)

lines_ffdihedral = 14233

print("\n")
print("Dihedral Potential Parameters")
for it_1 in range(0, len(dihedrals)):
    ID1 = atomtype[int(dihedrals[it_1,2]-1)]
    ID2 = atomtype[int(dihedrals[it_1,3]-1)]
    ID3 = atomtype[int(dihedrals[it_1,4]-1)]
    ID4 = atomtype[int(dihedrals[it_1,5]-1)]
    # ----------------------------------------------------------------------------
    ofi = open("ffdihedral", 'r')
    for it_2 in range(0, lines_ffdihedral):
        dump = ofi.readline()
        count = 0
        flag = 0
        ID1_fileff = 0
        ID2_fileff = 0
        ID3_fileff = 0
        ID4_fileff = 0
        for it_3 in range(0,len(dump)):
            if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                flag = 0
            if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                # ------------------------------------------------------------
                if (flag == 0):
                    count = count + 1
                    value = str()
                    flag = 1
                # ------------------------------------------------------------
                value = value+''.join(dump[it_3])
                
                if it_3 != (len(dump)-1):
                    if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t'):
                        if count == 1:
                            ID1_fileff = value
                        if count == 2:
                            ID2_fileff = value
                        if count == 3:
                            ID3_fileff = value
                        if count == 4:
                            ID4_fileff = value
                            break
        if count == 4:
            if (ID1 == ID1_fileff) & (ID2 == ID2_fileff) & (ID3 == ID3_fileff) & (ID4 == ID4_fileff):
                bit1 = str()
                bit2 = str()
                bit3 = str()
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t') or (dump[it_4] == '\n'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 7) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 8) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                    if (count2 == 6) & (flag == 1):
                        bit3 = bit3+dump[it_4]
                # ---------------------------------------
                bit1 = float(bit1)
                bit2 = float(bit2)
                bit3 = float(bit3)
                bit1 = "{0:.5f}".format(bit1/4.184)
                bit2 = "{0:.5f}".format(bit2)
                bit3 = "{0:.5f}".format(bit3)
                bit1 = str(bit1)
                bit2 = str(bit2)
                bit3 = str(bit3)
                # ---------------------------------------
                out = str(int(dihedrals[it_1,1]))+str("\t")+bit1+str("\t")+bit2+str("\t")+bit3
                print(out)
            if (ID4 == ID1_fileff) & (ID3 == ID2_fileff) & (ID2 == ID3_fileff) & (ID1 == ID4_fileff):
                bit1 = str()
                bit2 = str()
                bit3 = str()
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t') or (dump[it_4] == '\n'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 7) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 8) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                    if (count2 == 6) & (flag == 1):
                        bit3 = bit3+dump[it_4]
                # ---------------------------------------
                bit1 = float(bit1)
                bit2 = float(bit2)
                bit3 = float(bit3)
                bit1 = "{0:.5f}".format(bit1/4.184)
                bit2 = "{0:.5f}".format(bit2)
                bit3 = "{0:.5f}".format(bit3)
                bit1 = str(bit1)
                bit2 = str(bit2)
                bit3 = str(bit3)
                # ---------------------------------------
                out = str(int(dihedrals[it_1,1]))+str("\t")+bit1+str("\t")+bit2+str("\t")+bit3
                print(out)

# Same comment that I wrote in the end of the section where I print the bond
# potential parameters hold here for the dihedrals.

# This part is only necessary if parameters for a dihedral type is not found: not
# finding it would mean that either it doesnt exist listed on the file (this
# could be due to you messing up the declaration of atom names in the beggin-
# ing of the code) OR if the dihedral is formed by a wildcard, case in which 
# you should check the parameters manually in ffdihedral file for the given 
# sequence of "atom names". This should be used in synergy with the dihedral_
# solution_CGenFF.py code, which will tell you the index of the dihedrals for
# which no potential exists explicitly for the given 4 atom types forming it
# in the ffdihedral file.
print("\n")
print("Dihedral atom names")
for it_1 in range(0, len(dihedrals)):
    # These receive the line indexes of the information that concerns the dihedraled a-
    # toms in the atomtype array.
    ID1 = atomtype[int(dihedrals[it_1,2]-1)]
    ID2 = atomtype[int(dihedrals[it_1,3]-1)]
    ID3 = atomtype[int(dihedrals[it_1,4]-1)]
    ID4 = atomtype[int(dihedrals[it_1,5]-1)]
    print(it_1 +1, ID1, ID2, ID3, ID4)
    
# ----------------------------------------------------------------------
# Unlike for bonds, angles and dihedrals, for improper I will simply read the
# due LigParGen file section and output the CGenFF atom types forming the im-
# proper below to then find the parameters in the bottom of the ffbonded.itp
# file. Remember that as before you may need to take into account impropers
# that are present only in the dra_bonded.itp file.
impropers = []
ofi = open("impropers", 'r')
for it_1 in range(0, number_of_impropers):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(6)):
            impropers.append(float(e))
impropers = np.array(impropers,float)
impropers = impropers.reshape(number_of_impropers,6)

print("\n")
print("Improper atom names")
for it_1 in range(0, len(impropers)):
    ID1 = atomtype[int(impropers[it_1,2]-1)]
    ID2 = atomtype[int(impropers[it_1,3]-1)]
    ID3 = atomtype[int(impropers[it_1,4]-1)]
    ID4 = atomtype[int(impropers[it_1,5]-1)]
    print(ID1, ID2, ID3, ID4)

# -----------------------------------------------------------------------------------
# Within the scope of the methodology that I use to find the parameters for each bond,
# angle and dihedral type, I am assuming each bond, angle and dihedral declared in the
# respective sections have a unique type (which in this case would match the ID). This
# seems to be the case for many of the compounds I was initially testing LigPargen for,
# and so I ultimately assumed it is a universal conclusion (i.e. would hold for all mo-
# lecules). If this happens, for sure I wont have any problem in the methodology implied 
# in the code to find the parameters. If not, there may be problems. I am going to make 
# therefore a sanity check on this here and print a warning in case more than one ID 
# shares a type.
for it_1 in range(0, len(bonds)):
    if bonds[it_1,0] != bonds[it_1,1]:
        print("BOND TYPE DOES NOT MATCH BOND ID")
        
for it_1 in range(0, len(angles)):
    if angles[it_1,0] != angles[it_1,1]:
        print("ANGLE TYPE DOES NOT MATCH ANGLE ID")
        
for it_1 in range(0, len(dihedrals)):
    if dihedrals[it_1,0] != dihedrals[it_1,1]:
        print("DIHEDRAL TYPE DOES NOT MATCH DIHEDRAL ID")