% function [ result, f] = membrane_current_Kv(params, rho, V)
%
% 	Current flow across membrane for GHK flux with Boltzmann open probability
%
% 	Inputs
%		params -> state of system
%                   params.membrane.V_mid = voltage in V at which open probability po = 50%
%                   params.membrane.V_width = voltage range over which open
%                   probability changes by ~23% (i.e. 27% or 73% rather than 50%)
%                   i.e. po(Vm) = 1 ./ (1 + exp( (params.membrane.V_mid-Vm)/params.membrane.V_width ) ); 
%
%		rho -> 	  optional density if different from params.rho (mM = moles/m^3)
%		V ->	  optional voltages if different from params.rho (V)
%
% 	Output
%		result -> Molecular current carried by each species (moles/sec)
%       f -> flux of water crossing membrane (m^3/sec)

function [ result, f] = membrane_current_Kv(params, rho, V)

	if (nargin<2) % Concentrations in each compartment
		rho = params.rho;
	end
	if (nargin<3) % Voltage in each compartment
		V = params.V;
	end

	Nm = params.membrane.N; % Membrane is between compartments n = Nm and n = Nm+1
	rho_a = rho(:,Nm); 	% Concentration on Side A
	rho_b = rho(:,Nm + 1);	% Concentration on Side B 
	Vm = V(Nm) - V(Nm+1) ; % Voltage from A to B
	Am = params.A(Nm+1); % Area of membrane

	po = 1 ./ (1 + exp( (params.membrane.V_mid-Vm)/params.membrane.V_width ) ); % Open probability

	u = params.z*params.F*Vm/(params.R*params.T) ; 
	uz = (abs(u)<1e-6); % Uncharged particles and zero membrane voltge

	result = po * Am * params.membrane.perm .* (u + uz) ./ (1 - exp(-u) + uz) .* (rho_a - rho_b.*exp(-u)); % Molecular current carried by each species (moles/sec)    
    f = params.membrane.Pf * Am * params.Vwater * ( sum(rho_b) - sum(rho_a)  ); %       flux of water crossing membrane (m^3/sec)
    