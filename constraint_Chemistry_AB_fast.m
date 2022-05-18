%function result = constraint_Chemistry_AB_fast(params, params_i, delta_t)
%
% 	Constraint on chemical reactions for a solution containing
%	a simple buffer, A + B <-> AB with fast rate and equilibrium
%              A * B = params.Kd * AB
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

function result = constraint_Chemistry_AB_fast(params, params_i, delta_t)
    
	rr = params.rr; 	% Reaction Rates
	rho = params.rho;	% Concentrations
    result =  -1 * rr;  % Default is reaction rate should be zero.
    tscale = 1; 		% Timescale in seconds
    
	% Now evaluate buffer reactions
	idx_A =  1; idx_B =  2; idx_AB = 3; 
    
    result(idx_A,:)  = ( max(0, rho(idx_AB,:)) - max(0, rho(idx_A,:)) .* max( 0, rho(idx_B,:)) / params.Kd ) / tscale; % Rate that will be zero at equilibrium
    result(idx_B,:)  = rr(idx_B,:)  - rr(idx_A,:); % Creation rate of B must much that of A
    result(idx_AB,:) = rr(idx_AB,:) + rr(idx_A,:); % Creation rate of AB must be inverse of A
    