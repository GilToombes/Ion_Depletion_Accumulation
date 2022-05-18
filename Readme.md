# Ion_Accumulation_Depletion


Ion_Accumulation_Depletion is a small library of Octave/Matlab functions for calculating ion and voltage concentrations during patch clamp experiments.

Examples are in the Examples folder.

## Function List

calculate_PGHK_from_I.m
    Calculate the GHK membrane permeability to each ion
    corresponding to a given current of each ion.
    
calculate_tip_diameter.m
    Calculate pipette tip diameter for a given pipette resistance

constraint_Chemistry_AB.m
 	Constraint on chemical reactions for a simple buffer, A + B <-> AB 

constraint_Chemistry_AB_AC.m
 	Constraint on chemical reactions for solution with 2 simple buffers
    A + B <-> AB  and   A + C <-> AC 

constraint_Chemistry_AB_fast.m
 	Constraint on chemical reactions for a solution containing 
    a simple buffer  A + B <-> AB where the buffer reaction 
    rates are much faster than all other dynamics in system (i.e. fast buffer)

constraint_Chemistry_HB.m 
	Constraint function for a fast HB proton buffer
            H+ + B <-> HB
 
constraint_Chemistry_HOHB.m 
	Constraint function for a HB proton buffer
	which included OH-
          H+ + B   <-> HB ;
          H+ + OH  <-> H2O;

constraint_J.m
	Constraint function for individual molecular fluxes

constraint_Jsol.m
	Constraint function for solvent flow
    Assumes all comparments are incompressible

constraint_Lambda.m
    Constraint function for compartment volume
    assuming compartment is infinitely deformable

constraint_Rho.m
	Constraint function for molecular concentrations

constraint_V.m
	Constraint function for voltage

create_protocol.m
    Function to set up all necessary variables for an experimental protocol

evolve_backward_euler.m
	Evolve state using the backward Euler method

initialize_system.m
	Initialize state to default values.

membrane_current_constant.m
	Constant flux at membrane - can give unphysical answers.

membrane_current_GHK.m
	GHK flux for membrane

membrane_current_Kv.m
	GHK flux for membrane with voltage-dependent membrane permeability

plot2svg.m
	Save figures to SVG function from Juerg Schwizer 23-Oct-2005

resistance.m
	Calculate electrical resistance of configuration
	
save_figure.m
	Save and scale figures to multiple formats
