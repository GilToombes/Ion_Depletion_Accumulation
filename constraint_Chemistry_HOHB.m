%function result = constraint_Chemistry_HBOH(params, params_i, delta_t)
%
% 	Constraint on chemical reactions for a solution containing
%	a simple buffer, H+ + B- <-> HB with equilibrium constant Ka, params.K_HB
%	The creation of OH- is also considered, and the buffer is assumed to have
%	very fast association/disassociation kinetics
%		
%	Solution composition
%		1 = H+, 2 = B-, 3 =  HB, 4 = OH-, 5 ... N -> other components without chemistry
%
%	Inputs:
%		params ->  current state
% 		params_i ->initial state (not used)
%		delta_t ->  timestep
%
%	Output : Error in reaction rates relative to those needed to ensure that
%		  1. Creation Rates of H+ matches sum of creation rates for B- and OH- 
%		  2. Creation Rates for B- and HB are inverse of each other 
%		  3. [H+] * [B+]/Ka = [HB]
%		  4. [H+] * [OH-]/K_OH = 1;
%		  5. 0 for all non-reactive species
%		


function result = constraint_Chemistry_HBOH(params, params_i, delta_t)
    
	rr = params.rr; 	% Reaction Rates
	rho = params.rho;	% Concentrations
        tscale = 1; % Timescale in seconds

    	result = -1 * rr; 	% Default is reaction rate should be zero.
    
	% Now evaluate buffering reactions
	idx_H =  1; idx_B =  2; idx_HB = 3; idx_OH = 4;
	result(idx_H,:) = rr(idx_H,:) - rr(idx_B,:) - rr(idx_OH,:); % H+ created with OH- and B-
        result(idx_B,:) = rr(idx_B,:) + rr(idx_HB,:); % Interconversion of H+ (and B-) <-> to HB
    	result(idx_HB,:) = (rho(idx_H,:) .* rho(idx_B,:) / params.K_HB  - rho(idx_HB,:)) / tscale; %  [H+] * [B-] / K_HB = [HB]
    	result(idx_OH,:) = (rho(idx_H,:) .* rho(idx_OH,:)/ params.K_HOH -1)/tscale; % Equilibrium of H+ and OH- with H20
    
    
    