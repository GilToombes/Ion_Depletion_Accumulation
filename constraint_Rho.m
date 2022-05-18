% function result = constraint_Rho(params_f, params_i, delta_t)
%
% 	Constraint function for molecular concentrations, rho
%
%	Inputs -
%		params_f -> current system state
%		params_i -> initial system state 
%		delta_t  -> time step between initial and current state
%
%	Outputs - 
%		result -> params_f.rho - ( params_i.rho * Lambda_i/Lambda_F + drho_dt * delta_t)
%

function result = constraint_Rho(params_f, params_i, delta_t)

    % Constraint on each density value.
    
    Nseg = params_f.Nseg; 	% Number of compartments.
    J = params_f.J; 		% Rate of molecules flowing between compartments
    Nrho = size(J,1); 		% Number of molecular components
    
    % Flux entering each compartment.
    flux = [zeros(Nrho,1), J(:,[1:(Nseg-2)]) - J(:,[2:(Nseg-1)]), zeros(Nrho,1)]; 
        % Since the flux is only specified for one boundary of each end compartment,
        % we assume zero flux, which would keep the concentration constant.
    
    % Volume
    Lambda_f = ones(Nrho,1) * params_f.Lambda; % Final volume of each compartment.
    Lambda_i = ones(Nrho,1) * params_i.Lambda; % Initial volume of each compartment.
    
    % Constraint on each density
    if (delta_t ~= inf)   
		% Finite time step
        	result = params_f.rho - (Lambda_i .* params_i.rho./Lambda_f + (flux./Lambda_f + params_f.rr) * delta_t);% Error in density
    else
        	% Infinite time step
	        delta_t = 0.1; % Calculate concentration change due to flux only over this interval.
        	result = -1 * (flux./Lambda_f + params_f.rr) * delta_t;% Error in density
    end
    
    % Correct for ends.
    result(:,params_f.n_bath) =  params_f.rho(:,params_f.n_bath) - params_f.rho_bath; % Concentration in bath
    result(:,params_f.n_pip) =   params_f.rho(:,params_f.n_pip) - params_f.rho_pip;  % Concentration in pipette

