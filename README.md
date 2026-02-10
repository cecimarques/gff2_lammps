(1) Make sure also you put the code for the dihedrals and duly talk about the need to tune more than one potential over the same dihedrals. Make sure also you say the code needs to be interfaced w/ the input file generated w/ the LigParGen (like I did for gff_lammps)

SM

(1) A sample of LAMMPS input scripts for OPLS-AA and CGenFF as well as microstates collected during the overall simulation are made available amongst the files of this SM
- put one lammps.in for CGenFF and one for OPLS-AA (pressing part)
- create a folder with *all* final_pos.dat files for each compound within both force fields so that the person can access the parameters ! PS: this should match the simulations in table S1.
- it might be worth also putting some of the step_xxx files !

(2) Average vlues of molar density, varsigma and upsilon (and standard error) for all compounds at all \(P_{zz}\) values considered in the MD simulation is made available as further supplemental material. The result of the FFT of the VACF at specific pressure conditions and the hydrostatic limits (average value only, to avoid confusion w/ ppl thinking we ran several simulations) for all compounds studied in our work is also made available. 
