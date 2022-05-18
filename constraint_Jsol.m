%function result = constraint_Jsol(params, params_i, delta_t)
%
% 	Constraint function for solvent flow, Jsol
%		This version assumes all compartments are incompressible (i.e. Jsol should be the same for all segments)
%		Assumes flow across membrane can be calculated with ideal solution equation.
%		
%	Inputs -
%		params   -> current system state
%		params_i -> initial system state -> not actually used
%		delta_t  -> time step between initial and current state
%
%	Outputs - 
%		result(j) -> Jsol(j)   - Jsol(j-1) for compartments to the right of a membrane
%			         Jsol(j+1) - Jsol(j)   for compartments to the left of a membrane
%			         Jsol(j)- Jsol_membrane for compartments next to membrane
%


function result = constraint_Jsol(params, params_i, delta_t)
	
	% Solvent flow at membrane 
	Nm = params.membrane_N; 	% Position of membrane
	rho_a = params.rho(:,Nm); 	% Concentration on Side A
	rho_b = params.rho(:,Nm + 1);	% Concentration on Side B 
	Am = params.A(Nm+1);            % Area of membrane
        
    Jsol_memb = params.membrane_Pf * Am * params.Vwater * ( sum(rho_b) - sum(rho_a)  ); % Solvent flow across membrane

	Np = params.Nseg; 		% Number of compartments
    Jsol = params.Jsol;		% Solvent flux
    
    % Fluxes for compartments to right of membrane
	if (Nm < (Np-1) ) 
         	idxr = [(Nm+1):(Np-1)];    
	        result(idxr) = Jsol(idxr) - Jsol(idxr-1);
    end
    
    % Flux for compartments on each side of membrane
	result(Nm) = Jsol(Nm) - Jsol_memb;
    
    % Flux for compartments on left side of membrane
	if (Nm > 1)
        	idxl = [1:(Nm-1)];
	        result(idxl) = Jsol(idxl+1) - Jsol(idxl);
	end
    
    
   
        
    
    