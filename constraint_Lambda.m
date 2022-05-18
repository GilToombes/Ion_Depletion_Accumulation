%function result = constraint_Lambda(params, params_i, delta_t)
%
% 	Constraint function for volume changes due to solvent flow
%       Assumes compartments are completely deformable.
%       Should only be applied to a completely deformable cell compartment
%		
%	Inputs -
%		params   -> current system state
%		params_i -> initial system state -> not actually used
%		delta_t  -> time step between initial and current state
%
%	Outputs - 
%       result -> 
%               result(j) = 0 for   j = 1, params_f.Nseg (bath and pipette electrode compartment volumes do not change
%               all other j               
%               params.Lambda(j) - ( params_i.Lambda(j) + (params.Jsol(j-1)- params.Jsol(j)) * delta_t)
%
%   Note if delta_t is infinite then a completely deformable cell can
%   change by an arbitrary amount.  To allow delta_t = inf to work for 
%   equilibrating concentrations must pin volume.
% 

function result = constraint_Lambda(params, params_i, delta_t)

    result = zeros(1, params.Nseg);   
	idx = [2:(params.Nseg-1)];
    
    % Constraint on each compartment volume
    if (delta_t ~= inf)   
		% Finite time step
        	result(idx) = params.Lambda(idx) - (params_i.Lambda(idx) + (params.Jsol(idx-1)- params.Jsol(idx)) * delta_t); % Error in volume
    else
        % Lock volume for infinite time step
	        result(idx) = params.Lambda(idx) - params_i.Lambda(idx);
    end
    
    
    
   
        
    
    