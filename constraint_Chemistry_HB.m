%function result = constraint_Chemistry_HB(params, params_i, delta_t)
%
% 	Constraint on chemical reactions for a solution containing
%	a simple buffer, H+ + B- <-> HB with Ka, params.Ka
%	which is assumed to have very fast association/disassociation kinetics
%
%	Solution composition
%		1 = H+, 2 = B-, 3 = HB, 4 ... N -> other components without chemistry
%
%	Inputs:
%		params ->  current state
% 		params_i ->initial state (not used)
%		delta_t ->  timestep
%
%	Output : Error in reaction rates relative to those needed to ensure that
%		  1. Creation Rates of H+ and B- are equal
%		  2. Creation Rates for B- and HB are inverse of each other 
%		  3. [H+] * [B+]/Ka = [HB]
%		  4. 0 for all non-reactive species
%		

function result = constraint_Chemistry_HB(params, params_i, delta_t)
    
	rr = params.rr; 	% Reaction Rates
	rho = params.rho;	% Concentrations
    	tscale = 1; 		% Timescale in seconds
    
    	result = -1 * rr;       % Default is reaction rate should be zero.
    
	% Now evaluate buffer reactions
	idx_H =  1; idx_B =  2; idx_HB = 3; 
	result(idx_H,:) = rr(idx_H,:) - rr(idx_B,:); % H+ created with B-
    result(idx_B,:) = rr(idx_B,:) + rr(idx_HB,:); % Interconversion of H+ (and B-) <-> to HB
    result(idx_HB,:) = (rho(idx_H,:) .* rho(idx_B,:) / params.K_HB  - rho(idx_HB,:)) / tscale; %  [H+] * [B-] / K_HB = [HB]

    
    
    