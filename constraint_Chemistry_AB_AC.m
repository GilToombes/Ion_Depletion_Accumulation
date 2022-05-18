%function result = constraint_Chemistry_AB_AC(params, params_i, delta_t)
%
% 	Constraint on chemical reactions for a solution containing up to 2
% 	buffers of species A. 
%                           A + B <-> AB with rates
%             forward reaction rate = params.Kon_AB * A * B
%             reverse reaction rate = params.Koff_AB * AB
%                           A + C <-> AC with rates
%             forward reaction rate = params.Kon_AC * A * C
%             reverse reaction rate = params.Koff_AC * AC
%
%	Solution composition
%		1 = A, 2 = B, 3 = AB, 4 = C, 5 = AC, 6 ... N -> other components without chemistry
%
%	Inputs:
%		params ->  current state
% 		params_i ->initial state (not used)
%		delta_t ->  timestep
%
%	Output : Error in reaction rates relative to those needed to ensure that
%		  1. rr(A)  - Creation Rate for A
%         2. rr(B)  - Creation Rate for B
%         3. rr(AB) - Creation Rate for AB
%		  4. rr(C)  - Creation Rate for C
%         5. rr(AC) - Creation Rate for AC
%	         rr - 0  for all non-reactive species

function result = constraint_Chemistry_AB_AC(params, params_i, delta_t)
    
	rr = params.rr; 	% Reaction Rates
	rho = params.rho;	% Concentrations
    result =  rr;       % Default is reaction rate should be zero.
    
	% Now evaluate buffer reactions
	idx_A =  1; idx_B =  2; idx_AB = 3; idx_C =  4; idx_AC = 5; 
    
    crate_AB = params.Koff_AB * max(rho(idx_AB,:), 0) - params.Kon_AB * max(rho(idx_A,:),0) .* max(rho(idx_B,:), 0) ;
    crate_AC = params.Koff_AC * max(rho(idx_AC,:), 0) - params.Kon_AC * max(rho(idx_A,:),0) .* max(rho(idx_C,:), 0) ;
    result(idx_A,:)  = rr(idx_A,:)  - crate_AB - crate_AC; 
    result(idx_B,:)  = rr(idx_B,:)  - crate_AB;
    result(idx_AB,:) = rr(idx_AB,:) + crate_AB;
    result(idx_C,:)  = rr(idx_C,:)  - crate_AC;
    result(idx_AC,:) = rr(idx_AC,:) + crate_AC;
    
    
    
    