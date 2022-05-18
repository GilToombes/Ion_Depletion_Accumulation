%function result = constraint_Hv1_State(params, params_i, delta_t)
%
% 	Constraint on state of Hv1
%
%	Inputs:
%		params ->  current state
% 		params_i ->initial state (not used)
%		delta_t ->  timestep
%
%       Key fields
%           params.membrane_Hv1_state -> column vector with occupancies of each
%                                        channel state
%           params.membrane_Hv1_alpha   -> transition rates at 0mV
%           params.membrane_hv1_z       -> charge with transition
%
%	Output : 
%         result = (si + rate * delta_t) - sf;
%

function result = constraint_Hv1_State(params, params_i, delta_t)
    
    r = params.membrane_Hv1_rate;  % Transition Rates
    z = params.membrane_Hv1_z;    % Charge of Transition
    Nm = params.membrane_N;	% Position of membrane
	Vm = params.V(Nm) - params.V(Nm+1) ; % Voltage from A to B

    sf = params.membrane_Hv1_state;     % Final State
    si = params_i.membrane_Hv1_state;   % Initial State

    % Intra-cellular proton concentration
    gamma = params.membrane_Hv1_gamma ;
    idx_H = 1; % Index for protons
    rhoH_pH6p5 = 10^(-6.5+3); % Concentration of protons at pH 6.5
    rp = params.rho(1,Nm) / rhoH_pH6p5; % proton concentration relative to pH 6.5
    
    P = exp( z * Vm * params.F / (params.R * params.T)).* r .* rp.^gamma ; % Calculate off-diagonal rates
    P = P - diag(sum(P,2));
    
    
    
    if isinf(delta_t)
        tau = 1; % Timescale of 1 second
        result = sf * P * tau;
    else
        result = si + sf * P * delta_t - sf;
    end
    result(1) = result(1) + (1-sum(sf));
    
    
    

    
    
    
    
    