% function result = membrane_current_GHK(params, rho, V)
%
% 	Molecular flow across membrane for GHK flux equation.
%
% 	Inputs
%		params -> state of system
%		rho -> 	  optional density if different from params.rho (mM = moles/m^3)
%		V ->	  optional voltages if different from params.rho (V)
%
% 	Output
%		result -> Molecular current carried by each species (mole/sec)

function result = membrane_current_GHK(params, rho, V)

	if (nargin<2) % Concentrations in each compartment
		rho = params.rho;
    end
    
	if (nargin<3) % Voltage in each compartment
		V = params.V;
	end

	Nm = params.membrane_N;	% Membrane is between compartments n = Nm and n = Nm+1
	rho_a = rho(:,Nm); 	% Concentration on Side A
	rho_b = rho(:,Nm + 1);	% Concentration on Side B 
	Vm = V(Nm) - V(Nm+1) ; % Voltage from A to B
	Am = params.A(Nm+1);  % Membrane Area

	u = params.z*params.F*Vm/(params.R*params.T) ; 
	uz = (abs(u)<1e-6); % Uncharged particles and zero membrane voltge

	result = Am * params.membrane_perm .* (u + uz) ./ (1 - exp(-u) + uz) .* (rho_a - rho_b.*exp(-u)); % Molecular current carried by each species (mole/sec)
