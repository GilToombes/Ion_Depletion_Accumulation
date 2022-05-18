%function result = constraint_V(params, params_i, delta_t)
%
% 	Constraint function for Voltage
%
%	Inputs -
%		params   -> current system state
%		params_i -> initial system state 
%		delta_t  -> time step between initial and current state
%
%	Outputs - 
%		result -> Charge - Capacitance * Voltage 
%
     

function result = constraint_V(params, params_i, delta_t)
    
	Np = params.Nseg; 	% Number of compartments
	Cl = params.Cl; 	% Longitudinal capacitance
    Ct = params.Ct; 	% Transverse capacitance
    V = params.V; 		% Voltage
    Lambda = params.Lambda; % Volume of compartments
	Q = sum((params.F * params.z * Lambda ) .* params.rho , 1);    	% Charge per compartment

	% Capacitative Charge per compartment
	idx = [1:(Np-1)];
    CV = [Cl(idx).*(V(idx)-V(idx+1)), 0] + [0, Cl(idx).*(V(1+idx)- V(idx))] + Ct.*V;

    % Scale factor for each element, excluding ends.
    sf = max( params.F * Lambda , [Ct(1), Cl(idx) + Cl(idx+1) + Ct(idx+1)]);
    
    % Calculate error
	result = (Q - CV)./sf;

	% Correct for end compartments which contain electrodes
	if ~isfield(params,'Rseries')
             % No series resistance compensation
        	 V_com = params.V_hold * (1 - 2*(params.n_pip>1)); % Command voltage at pipette without series resistance compenstaion
    else
            % Series Resistance Compensation
         	Npip_current = max(params.n_pip-1,1); % Index of current fluxes at pipette end     
            I_pip = params.F * sum(params.J(:,Npip_current) .* params.z); % Current flowing at pipette electrode
         	V_com_i = params_i.V(params.n_pip); % Previous command volutage.
            if (nargin>2) 
                    % Calculate percent change of command voltage during time step of delta_t
                    % params.Rs_tau gives the time constant for the series resistance compensation    
                fchange = 1 - params.Rseries_tau / (params.Rseries_tau + delta_t);
            else
                fchange = 1;
            end
         	V_com = V_com_i * (1-fchange) + fchange * (params.V_hold + I_pip * params.Rseries) * (1 - 2*(params.n_pip>1)); % Command voltage at pipette without series resistance compensation
	end
	result(params.n_bath) = (V(params.n_bath) - 0); % Bath voltage clamped at 0
    result(params.n_pip) = (V(params.n_pip) - V_com); % Pipette voltage clamped at command voltage 


        