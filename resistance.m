% function result = resistance(params)
%
% 	Calculate resistance 
%
% 	Inputs -> params - describes system
%
% 	Output -
% 		If no membrane 
%			result = resistance between left and right electrodes
%	 	If membrane
%	 		result(1) = resistance between left electrode and membrane in Ohms
%			result(2) = resistance between membrane to right electrode in Ohms

function result = resistance(params)

	Np = params.Nseg;
	idx = [1:(Np-1)];
	
	% Interface Area
	Ai = sqrt(params.A(idx) .* params.A(idx+1) ); 

	% Width	
	dx = params.x(idx+1) - params.x(idx);

	% Drift Current
	sigma = sum( (params.rho(:,idx+1)+params.rho(:,idx))/2  .* (  (params.D .* params.z .* params.z) * ones(1,Np-1) ) , 1) * params.F^2 / (params.R * params.T);

	dR = dx ./ (Ai .* sigma);
    
	if isfield(params,'membrane_func')
		Nm = params.membrane_N;
		if (Nm>1)
			result(1) = sum( dR([1:(Nm-1)]) );
		else
			result(1) = 0;
		end

		if (Nm<(Np-1))
			result(2) = sum(dR([(Nm+1):(Np-1)]) );
		else
			result(2) = 0;
		end
	else
		result = sum(dR);
	end