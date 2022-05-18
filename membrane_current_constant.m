% function result = membrane_current_constant(params, rho, V)
%
% 	Molecular flow across membrane for constant current.
%   Note need to be cautious using this as large currents can cause
%   unphysical (i.e. negative) concentrations.
%
% 	Inputs
%		params -> state of system
%                 params.membrane_J = [ J_a, J_b, J_c, ...]';  
%                               Molecular current carried by each species
%		rho -> 	  optional density if different from params.rho
%		V ->	  optional voltages if different from params.rho
%
% 	Output
%		result -> Molecular current carried by each species

function result = membrane_current_constant(params, rho, V)

	result = params.membrane_J ; 