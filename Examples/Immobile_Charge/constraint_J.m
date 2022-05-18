%function result = constraint_J(params, params_i, delta_t)
%
% 	Constraint function for molecular fluxes, J
%
%	Inputs -
%		params   -> current system state
%		params_i -> initial system state -> not actually used
%		delta_t  -> time step between initial and current state
%
%	Outputs - 
%		result -> J - params.J
%			  J 		the calculated molecular flux due to diffusion, drift and convection
%			  params.J 	the current molecular flux (which should be equal to J).			
%

function result = constraint_J(params, params_i, delta_t)

    Nseg = params.Nseg;      % Number of segments
    Ncom = length(params.z); % Number of components
    rho = params.rho;        % Concentration of each component
    V = params.V; 	     % Voltage
    Jsol = params.Jsol;      % Solvent Flow Rates
    
    idx = [1:(Nseg-1)];	     			  % Index for calculating differences
    dx = params.x(idx+1) - params.x(idx); 	  % Distance between each pair of segments
    Ai = sqrt(params.A(idx) .* params.A(idx+1) ); % Area at interface between each pair of segments
    rho_avg = (rho(:,idx+1)+rho(:,idx))/2; 	  % Density at interface
    T = params.D * (Ai./dx ); 			  % Diffusive transport coefficient

    % Diffusive flux
    J_diff = -1 * T .* (rho(:,idx+1) - rho(:,idx));
    
    % Drift flux
    sf = -1 * params.F / (params.R * params.T) ; 
    J_drift =  sf * (params.z * (V(idx+1)-V(idx))) .* T .* rho_avg;

    % Convective Flux
    J_conv = ((params.D > 0) * Jsol) .*  rho_avg ;     
	% Allow for non-mobile ions by giving them a diffusion coefficient of 0.

    % Total Flux
    J = J_diff + J_drift + J_conv ;
    
    % Flux for segments separated by membrane
    if isfield(params,'membrane_func')
        J(:,params.membrane_N)= params.membrane_func(params);
    end
       
    % Flux for segments separated by patch
    if isfield(params,'patch_func')
        J(:,params.patch_N) = params.patch_func(params);
    end
    
        
    % Calculate difference between calculated and current flux 
    result = J - params.J;
    
