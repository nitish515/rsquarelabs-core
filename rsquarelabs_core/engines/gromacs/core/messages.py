from settings import version, weblink

welcome_message = """
********************************************
[  RSQUARE-LABS GROMACS AUTOMATION SCRIPT  ]
********************************************
Version """ + version + """
Check out """ + weblink + """ for more information.

  ****CHEERS & HAPPY RESEARCH****

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License version 2 (GNU GPL v2)
Do not remove the original authors name & credits in redistributions or
modified version.
        """
backup_folder_already_exists = """
    POSSIBLE SOLUTION: A backup folder already present.
    You may loose your work if we replace , so change your
    Working Directory or delete the folder BACKUP
"""

# used by grompp in Protein Min Neutralise step
ions_mdp ="""
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nsteps		= 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
coulombtype	    = PME		; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; Short-range electrostatic cut-off
rvdw		    = 1.0		; Short-range Van der Waals cut-off
pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)
"""

# used by grompp in Protein Min Production step
minim_mdp = """
; minim.mdp - used as input into grompp to generate em.tpr
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nsteps		= 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
coulombtype	    = PME		; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; Short-range electrostatic cut-off
rvdw		    = 1.0		; Short-range Van der Waals cut-off
pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)
"""


write_em_mpd_data = """
; LINES STARTING WITH ';' ARE COMMENTS
title         = Minimization    ; Title of run

; Parameters describing what to do, when to stop and what to save
integrator    = steep        ; Algorithm (steep = steepest descent minimization)
emtol         = 100.0        ; Stop minimization when the maximum force < 1.0
                               kJ/mol
emstep        = 0.01         ; Energy step size
nsteps        = 50000        ; Maximum number of (minimization) steps to
perform
energygrps    = system       ; Which energy group(s) to write to disk

; Parameters describing how to find the neighbors of each atom and how to
  calculate the interactions
nstlist       = 1            ; Frequency to update the neighbor list and long
                               range forces
ns_type       = grid         ; Method to determine neighbor list (simple, grid)
rlist         = 1.0          ; Cut-off for making neighbor list
                               (short range forces)
coulombtype   = PME          ; Treatment of long range electrostatic
                               interactions
rcoulomb      = 1.0          ; long range electrostatic cut-off
rvdw          = 1.0          ; long range Van der Waals cut-off
pbc           = xyz          ; Periodic Boundary Conditions (yes/no)
"""

create_em_mdp_data = """
; LINES STARTING WITH ';' ARE COMMENTS
title        = Minimization  ; Title of run

; Parameters describing what to do, when to stop and what to save
integrator     = steep       ; Algorithm (steep = steepest descent minimization)
emtol          = 100.0       ; Stop minimization when the maximum force < 1.0
                               kJ/mol
emstep         = 0.01        ; Energy step size
nsteps         = 50000       ; Maximum number of (minimization) steps to perform
energygrps     = Protein UNK ; Which energy group(s) to write to disk

; Parameters describing how to find the neighbors of each atom and how to
calculate the interactions
nstlist        = 1           ; Frequency to update the neighbor list and long
                               range forces
ns_type        = grid        ; Method to determine neighbor list (simple, grid)
rlist          = 1.0         ; Cut-off for making neighbor list (short range
                               forces)
coulombtype    = PME         ; Treatment of long range electrostatic
                               interactions
rcoulomb       = 1.0         ; long range electrostatic cut-off
rvdw           = 1.0         ; long range Van der Waals cut-off
pbc            = xyz         ; Periodic Boundary Conditions (yes/no)
"""