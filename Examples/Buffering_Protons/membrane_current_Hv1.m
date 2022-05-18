% function membrane_current_Hv1(params,rho,V)
%
% 	Molecular flow across membrane for model of Hv1
%
% 	Inputs
%		params -> state of system
%		rho -> 	  optional density if different from params.rho
%		V ->	  optional voltages if different from params.rho
%               
%
% 	Output
%		result -> Molecular current carried by each species

function result = membrane_current_Hv1(params, rho, V)

	if (nargin<2)
		rho = params.rho;
	end
	if (nargin<3)
		V = params.V;
	end

	Nm = params.membrane_N;	% Position of membrane
	rho_a = rho(:,Nm); 	% Concentration on Side A
	rho_b = rho(:,Nm + 1);	% Concentration on Side B 
	Vm = V(Nm) - V(Nm+1) ; % Voltage from A to B
	
    perm = params.membrane_perm ; % General membrane permeability
    % H+ current due to Hv1
    idx_H = 1; % H+ index
    idx_open = 4; % Open state of Hv1
    perm(idx_H) = perm(idx_H) + params.membrane_Hv1_state(idx_open) * params.membrane_Hv1_permH;

    Am = params.A(Nm+1);
    
	u = params.z*params.F*Vm/(params.R*params.T) ; 
	uz = (abs(u)<1e-6); % Uncharged particles and zero membrane voltge

	result = Am * perm .* (u + uz) ./ (1 - exp(-u) + uz) .* (rho_a - rho_b.*exp(-u));

