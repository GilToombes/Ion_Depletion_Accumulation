%function result = constraint_Chemistry_AB(params, params_i, delta_t)
%
% 	Constraint on chemical reactions for a solution containing
%	a simple buffer, A + B <-> AB with rates
%             forward reaction rate = params.Kon * A * B
%             reverse reaction rate = params.Koff * AB
%
%	Solution composition
%		1 = A, 2 = B, 3 = AB, 4 ... N -> other components without chemistry
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
%		  4. rr - 0  for all non-reactive species

function result = constraint_Chemistry_AB(params, params_i, delta_t)
    
	rr = params.rr; 	% Reaction Rates
	rho = params.rho;	% Concentrations
    result =  rr;       % Default is reaction rate should be zero.
    
	% Now evaluate buffer reactions
	idx_A =  1; idx_B =  2; idx_AB = 3; 
    
    crate = params.Koff * max(rho(idx_AB,:), 0) - params.Kon * max(rho(idx_A,:),0) .* max(rho(idx_B,:), 0) ;
    result(idx_A,:)  = rr(idx_A,:)  - crate; 
    result(idx_B,:)  = rr(idx_B,:)  - crate;
    result(idx_AB,:) = rr(idx_AB,:) + crate;
    
    
    
    