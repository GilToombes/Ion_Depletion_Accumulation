Simulation of proton channel currents flowing through voltage-gated proton channel, Hv1.
       Ideal : No ion accumulation or depletion
       Whole Cell :  Raccess = 5 MOhms;  Cell volume = 1 pL; Membrane Area = 1000 um^2
       MacroPatch :  Raccess = 0.23 Mohms; Membrane Area = 300 um^2
        pH 5.5 versus 7.5
        
Main Function
       Plot_Hv1_Giant_Patch_vs_wholecell_v2.m -> Does Simulation
       
Ancillary Functions
       constraint_Hv1_State -> characterizes voltage-dependent transitions of Hv1 channel
       membrane_current_Hv1 -> I-V relationship for Hv1 channel
       