% function params = initialize_system(params)
%
% 	Function for initializing all system parameters for a calculation
%
%   Input 
%           Structure, params, is optional as any non-specified fields are
%           replaced with defaults.
%
%           Setup          
%                 params.config  -> Can specify experimental configuration
%
%                        'outside-out' -> default
%                          Pipette electrode in first compartment (n=1)
%                          Pipette is divided into Nseg-1 compartments
%                          Bath electrode in final compartment   (n=Nseg)
%                       
%                        'whole-cell'
%                          Pipette elecrode in first compartment (n=1)
%                          Pipette is divided into Nseg-2 compartments
%                          Cell is compartment n = Nseg-1
%                          Bath is compartment n = Nseg
%       
%                        'inside-out'
%                          Bath electrode in first compartment (n=1)
%                          Pipette divided into Nseg-1 compartments
%                          Pipette electrode in last compartment n = Nseg
%
%          If using config = 'outside-out','whole-cell', or 'inside-out'
%
%              params.Nseg   ->  Number of compartments
%                                   Nseg = 11 is default
%
%              Pipette Shape (only useful for a conical pipette)
%                 params.theta  ->  Pipette tip angle in radians
%                                   theta = 10 * pi/180 
%                                        10 degrees by default
%                 params.r_tip  ->  Pipette tip radius in m 
%                                   r_tip = 0.5e-6 is default 
%                                       1 micron tip diameter by default
%                 params.r_pip  ->  Pipette inner radius at electrode in m
%                                   r_pip = 0.5e-3 is default
%                                       1mm inner diameter pipette glass 
%                 params.d_ratio -> % Inner pipette wall diameter / Outer Pipette wall Diameter
%                                   d_ratio = 0.75 is default
%
%               Solutions
%                 params.z    Nspecies x 1 vector of charge of each ionic species
%                               z = [1 -1]'; % Default is KCl
%                 params.D    Diffusion coefficient of each ionic species in m^2/s
%                               D = 1.957e-9* [1 1.038]';   Default is K+,  Cl-  
%                                    http://web.med.unsw.edu.au/phbsoft/mobility_listings.htm
%                                    http://pubs.acs.org/doi/pdf/10.1021/ja01622a016%
%                 params.rho_pip  : Nspecies x 1 vector of concentration of each ionic species in pipette solution (mM)
%                                 rho_pip = [150 150]'; 150mM KCl default
%                 params.rho_bath : Nspecies x 1 vector of concentration of each ionic species in bath solution (mM)
%                                 if not specified rho_bath = rho_pip;
%                                 Can also be a function handle to give concentration depending on time.
%                                 rho_bath = rho_bath_function(time)
%                                 then initialize to rho_b = params.rho_bath(params.t);
%
%               Voltages
%                   params.V_hold -> Holding voltage in Volts
%                                 V_hold = 0 ; default.
%                   params.Rseries -> Optional series resistance compensation in Ohms
%                   parmams.Rseries_tau -> Timescale for series resistance compensation (s)
%                               Rseries_tau = 1e-5; % Default is 10 microseconds
%                               % Stable when tau > Cpipette * Rpip                         
%
%               Membrane
%                params.membrane_func -> Optional function handle describing membrane
%                                      permeability
%                                    e.g. membrane_func = @membrane_current_GHK
%                                       for a GHK flux 
%                                      If no function is specified there is no membrane
%
%                params.membrane_x ->    Position of membrane in m
%                                        Not relevant for 'whole-cell'
%                                       e.g. membrane_x = 0 -> at pipette tip
%                                       e.g. membrane_x = -3e-6 -> 3 microns
%                                       from tip
%
%               Cell Parameters (for config = 'whole-cell')
%               params.A_cell  ->       Area of cell membrane in m^2
%                                       A_cell = 1e-9; % Default is 1000 um^2
%               params.Lambda_cell ->   Volume of cell in m^3 (
%                                  Lambda_cell = 1e-15 as default     (1pL)
%
%
%               params.Contraints         Structure dedscribing which parameters evolve
%                  Constraints.type = 'basic';              % Concentrations (rho), molecular fluxes (J), and voltage (V)
%                                   = 'basic_perfusion';    % Concentrations (rho), molecular fluxes (J), and voltage (V)
%                                                           % Allow concentrations of pipette and bath compartments to evolve
%                                   = 'solvent';            % Concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol)
%                                                            Constant compartment volumes - uniform water flow, Jsol, through all
%                                   = 'solvent_cell_volume' % Concentrations (rho), molecular fluxes (J), voltage (V), solvent flux at membrane (Jsol), and cell volume
%                                                           % Water flow across whole-cell membrane, but no water flow in pipette
%                                   = 'buffer';             % Concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered (rr)
%
%   Output 
%           Structure, params, containing fields.
%
%           Physical constants
%                 params.R      % Gas constant in J/(mol.K)
%                 params.F      % Faraday's constant in Coulombs per mole
%                 params.Vwater % Molar volume of water (in m^3/mol)
%                 params.z    Nspecies x 1 vector of charge of each ionic species
%                               z = [1 -1]'; % Default is KCl
%                 params.D    Diffusion coefficient of each ionic species in m^2/s
%                               D = 1.957e-9* [1 1.038]';   Default is K+,  Cl-  
%                                    http://web.med.unsw.edu.au/phbsoft/mobility_listings.htm
%                                    http://pubs.acs.org/doi/pdf/10.1021/ja01622a016%
%           State 
%                 params.t      % Current time in seconds
%                                   0s = Default
%                 params.T      % Temperature in Kelvin 
%                                   20C = 293K = Default
%
%           Amplifier
%                 params.V_hold -> Holding voltage in Volts
%                                 V_hold = 0 ; default.
%                 params.Rseries -> Optional series resistance compensation in Ohms
%                 parmams.Rseries_tau -> Timescale for series resistance compensation (s)
%                               Rseries_tau = 1e-5; % Default is 10 microseconds
%                               % Stable when tau > Cpipette * Rpip 
%
%           Physical Configuration
%                 params.Nseg   % Number of compartments
%                 params.n_pip  % Compartment containing pipette elecrode
%                                   n_pip = 1 for 'outside-out' and 'whole-cell'
%                                   n_pip = N_seg for 'inside-out'
%                 params.n_bath % Compartment containing bath electrode
%                                   n_bath = N_seg for 'outside-out' and 'whole-cell'
%                                   n_bath = 1 for 'inside-out'
%                 params.x      % Vector of positions (m) of the reference point within each compartment.
%                 params.A      % Vector of cross-sectional area (m^2) at the reference point within each compartment
%                 params.Lambda % Vecor of volumes (m^3) of each compartment
%                 params.rho    % Nspecies * Nseg matrix of concentrations (mM) of ionic species
%                                   params.rho(a,m) is concentration of
%                                   "a"-th specie in compartment m.
%                 params.J      % Nspecies * (Nseg-1) matrix of molecular fluxes (moles/sec)
%                                   params.J(a,m) is flow rate of species
%                                   "a" from compartment m into m+1
%                 params.J_sol  % 1 * (Nseg-1) matrix of solvent flow rate (m^3/sec)
%                                   params.J_sol(m) is flow rate of solvent from compartment m into m+1
%                 params.rr     % Nspecies * Nseg matrix of molecular creation rates (mM/sec)
%                                   params.rr(a,m) = rate species a in
%                                   compartment m is being created.
%
%           Membrane
%              params.membrane_func : Optional function handle describing membrane permeability
%                      e.g. membrane_func = @membrane_current_GHK for a GHK flux 
%                      If 'membrane_func' fiel,d is omitted there is no membrane
%              params.membrane_x ->   Optional position of membrane in m
%                                     Only needed if params.x/A/lambda are
%                                     not specified explicitly.
%                      e.g. membrane_x = 0 -> at pipette tip
%                      e.g. membrane_x = -3e-6 -> 3 microns from tip
%               params.membrane_N ->  Membrane is located between
%                           compartment n = membrane_N  and n = membrane_N+1
%
%           Capacitances
%               params.Cl(1 ... (Nseg-1) ) 
%                       Longitudinal capacitance describing charge 
%                       due to voltage differences between adjacent compartments.
%                       Charge in compartment m due to m/m+1 boundary 
%                       = params.Cl(m) * (params.V(m) - params.V(m+1)) 
%               
%                       For compartments separated by water 
%                               calculation assumes dielectric constant of 80 
%                       For compartments separated by membrane
%                               calculation assumes specific membrane capacitance of 
%                                   Cspecific = 1e-2 Farads/m^2
%
%               params.Ct(1 ... Nseg)
%                       Transverse capacitance describing charge
%                       due to voltage differences between compartments in pipette and bath.
%                       Charge in compartment m due at walls is
%                       = params.Ct(m) * (params.V(m) - params.V(params.Nbath))
%
%                       Calculation assumes pipette wall dielectric
%                       constant of 3
%  
%           Constraints -> Structure specifying constraint functions
%               Constraints.param     -> list of parameter names being constrained.
%                           e.g. param =      {'rho',        'J',           'V'};              basic
%                           e.g. param =	  {'rho', 		'J',            'V',	'Jsol',			'rr'};  buffer
%
%               Constraints.fhandle   -> function constraining each of listed parameters
%                           e.g. fhandle =    {@constraint_Rho, @constraint_J,          @constraint_V};  basic
%                           e.g. fhandle =    {@constraint_Rho, @constraint_J,          @constraint_V, @constraint_Jsol, 	@constraint_Chemistry_AB};       buffer
%
%               Constraints.segments  -> Compartments to which constraints are applied
%                           e.g. segments = { [2,params.Nseg-1],      [1,params.Nseg-1],      [1,params.Nseg],}; % basic
%                                             buffer and pipette compartment concentrations fixed 
%                                segments = { [2,params.Nseg-1],      [1,params.Nseg-1],      [1,params.Nseg],	[1,params.Nseg-1],	[2, params.Nseg-1]}; % buffered
%                                             buffer and pipette compartment concentrations fixed   
%    			
%               Constraints.components-> components of parameter to which constraints are applied
%                           e.g. components ={[1,size(params.rho,1)], [1,size(params.rho,1)], [1 1]}; % basic
%                           e.g. components ={[1,size(params.rho,1)], [1,size(params.rho,1)], [1 1],			[1 1],			[1 3]};  % buffered
%                                               Only first 3 components are buffered in AB buffer system
%
%               Constraints.stepsize  -> amount to vary parameter by when calculating derivatives
%                           e.g. stepsize = 	{0.001,                   0.01,                 1e-5};  % basic
%                           e.g. stepsize = 	{0.001,                   0.01,                 1e-5,			1e-18,			1};      % buffered
%                                               Concentrations -> df/drho = (f(rho + 0.001 mM, J, V, Jsol, rr) - f(rho, J, V, Jsol, rr)) / (0.001mM)
%                                               Flux              df/dJ =   (f(rho, J + 0.01 mole/sec,V,Jsol,rr) - f(rho,J,V,Jsol,rr)) / (0.01 mole/sec)                                                                
%                                               Voltage           df/dV =   (f(rho, J, V + 0.01mV, Jsol, rr) - f(rho, J, V, Jsol, rr)) / (0.01 mV/sec)
%                                               Jsol            df/dJsol =  (f(rho, J, V, Jsol + 1e-18 m^3/s, rr) - f(rho, J, V, Jsol, rr)) / (1e-18 m^3/s)                             
%                                               Reactions         df/drr =  (f(rho, J, V, Jsol, rr + 1 mole/sec) - f(rho, J, V, Jsol, rr)) / (1 mole/s)                             
%
%               Constraints.err_tol    -> residual error allowed after optimization.
%                           e.g. err_tol = 	{0.0000001,           1e-18,              1e-5}; % basic
%                                               0.1nM,    1attomole/sec,       10microVolts
%                           e.g. err_tol = 	{0.001,               1e-18,              1e-5,			1e-19,			1e-3};   % buffered
%                                               1microM,  1attomole/sec,      10microVolts,    0.1 fL/sec, 1 microMole/sec)
%
%               Constraints.idx_v    - Cell array of indices for each constraint function outputs
%                                     idx_v{u} list the constrained indices for constraint function Constraints.fhandle{u} 
%                                           i.e.         comparment is element of Constraints.components{u}, 
%                                                        component is element of Constraints.components{u}
%               Constraints.idx_c    - Cell array of indices mapping the constraint state vector to the constrained parameters
%                                     For the u-th constrained parameter, Constraints.param{u},
%                                     the n-the element, idx_c{u}(n), in contraint state vector maps to params.(Constraints.param{u})(idx_v{u}(n))  
%
%                                     note - The ordering of elements within idx_v{u}/idx_{c} is intentionally chose to minimize the bandwidth of the Hessian
    
function params = initialize_system(params)

    
    % Create parameters structure if none given.
    if (nargin<1) 
		params.R = 8.314;  
    end
    
    % Physical Constants
    if ~isfield(params,'R'); params.R = 8.314; end; % Gas Constant in J/(mole.K)
	if ~isfield(params,'F'); params.F = 96485; end; % Faraday's constant C/mole
	if ~isfield(params,'T'); params.T = 293 ;  end; % Temperature in Kelvin.
    if ~isfield(params,'Vwater'); params.Vwater = 18e-6; end % Molar volume of water (in m^3/mol)
    
    % Time
	if ~isfield(params,'t')
		params.t = 0;
    end	

    % Patch-Clamp Configuration 
	if ~isfield(params,'config'); 	
		params.config = 'outside-out';  % Default to outside-out
	end 
                                   
    % Pipette Shape
   	if (~isfield(params,'x'))|(~isfield(params,'A'))	
	    if ~isfield(params,'theta') ; params.theta = 10 * pi/180; end;  % Pipette tip angle in radians
        if ~isfield(params,'r_tip') ; params.r_tip = 0.5e-6; end; 	% Pipette Tip radius
        if ~isfield(params,'r_pip') ; params.r_pip = 0.5e-3; end; 	% Pipette inner radius at electrode
	end


    % Position of reference point in each compartment
    if ~isfield(params,'x')

		if ~isfield(params,'Nseg') params.Nseg = 11; end % Number of compartments

	        % Reciprocal distances 
		        
                % Old
                ua = sin( params.theta/2) / params.r_pip;
                ub = sin( params.theta/2) / params.r_tip;

                % Position reference points depending on configuration 

	        if ~isfield(params,'membrane_func');
        	    % No membrane
		            u = ua + [0:(params.Nseg-1)]/(params.Nseg-1) * (ub-ua);
			    params.x = 1/ub - 1./u;          
        	else
	            % Membrane
        	    switch  params.config

                	case {'whole-cell'}
    				% left side is pipette electrode, right side is bath electrode
        			% params.Nseg-1 is at x=0 is pipette tip
        			% params.Nseg is the bath -> locate at x = 1/ub;

	               		u = ua + [0:(params.Nseg-1)]/(params.Nseg-2) * (ub-ua);
        	        	params.x = 1/ub - 1./u; % Position of each 
                		params.x(params.Nseg) = 1/ub; % Bath position
                		params.membrane_N = params.Nseg-1; % compartment on left side of membrane
                
                	case {'outside-out'}
		                % left side is pipette electrode,
                        % x=0 on right side at pipette tip
                	    % membrane located some distance along the pipette							
    
                        %um = 1/ (params.membrane_x + 1/ub);
                        um = ub/ (params.membrane_x*ub + 1);
                        phi = (ub-um)/(ub-ua);
                        Nm =  floor(phi * (params.Nseg-1));			
                    	Nm2 = params.Nseg-(Nm+2);
                    	u =  [ (ua +[0:Nm2]/(Nm2 +(Nm2==0)) * (um -ua)), (um + ([0:Nm] + (Nm==0) )/(Nm + (Nm==0)) * (ub-um)) ];
                     	params.membrane_N = Nm2+1;
                    	params.x = 1/ub - 1./u;
                    
        	        case {'inside-out'}    
                		% x=0 on left side at pipette tip
                        % right side is pipette electrode
	                    % membrane located some distance along the pipette
		
                        um = 1/ (params.membrane_x + 1/ub);
		                phi = (ub-um)/(ub-ua);
                    	if (phi<0.5) 
                           		Nm =  ceil(phi * (params.Nseg-2));			
                    	else
                           		Nm =  floor(phi * (params.Nseg-2));
                    	end
                    	Nm2 = params.Nseg-(Nm+2);
                    	u =  [(ub + [0:Nm]/(Nm + (Nm==0)) * (um-ub)), (um +([0:Nm2] + (Nm2==0))/(Nm2 +(Nm2==0)) * (ua -um)) ];
                    	params.membrane_N = Nm+1;
                    	params.x = 1./u - 1/ub;
       				end              
        	end
    end

    
    	% Confirm number of elements matches
	    params.Nseg = length(params.x);
	    idx = [1:(params.Nseg-1)];	
        
	% Area at each point
    if ~isfield(params,'A')
    
            %Existing
            R = abs(params.x) * sin( params.theta/2) + params.r_tip; % Radius at each point
           	params.A = pi * R.^2 / cos(params.theta/4)^2; % Area at each point
            
            
        % Include cell Area if whole-cell
        if strcmp( params.config, 'whole-cell')
        	% R = (params.x-params.x(1)) / (params.x(params.Nseg-1) - params.x(1)) * (params.R_ends(2) - params.R_ends(1)) + params.R_ends(1); 		
			
            
            %params.A = pi * R.^2;
            params.A = pi * R.^2 / cos(params.theta/4)^2; % Area at each point
            
            % Cell Area
			if ~isfield(params,'A_cell')
				params.A_cell = 1e-9; % Default cell area in m^2 
			end
			params.A(params.Nseg) = params.A_cell;
        end
    end
            
	% Volume of each section
	   if ~isfield(params,'Lambda')	
		R = sqrt(params.A / pi);
		gamma = ( R(idx+1) - R(idx) ) ./ ( R(idx+1) + R(idx) ) ;
		Vf = [params.A(idx) .* (params.x(idx+1) - params.x(idx))/2 .* (1 - 2/3*gamma.^2 - gamma.^3 /3), 0];
		Vb = [0, params.A(idx+1) .* (params.x(idx+1) - params.x(idx))/2 .* (1 - 2/3*gamma.^2 + gamma.^3 /3)];

		% Cell Volume
		if strcmp(params.config,'whole-cell')
			if ~isfield(params,'Lambda_Cell')
				params.Lambda_Cell = 1e-15 ; % Default cell volume = 1fL
			end
			
			Vf(params.Nseg - 1 ) = params.Lambda_Cell;
			Vb(params.Nseg) = params.Lambda_Cell;
		end

		params.Lambda = Vf + Vb;
       end

	% Longitudinal capacitance
	   if ~isfield(params,'Cl');
		epsilon_sol = 80 * 8.85e-12; % Permittivity of Solution
		params.Cl = [ epsilon_sol * sqrt( params.A(idx) .* params.A(idx+1)) ./ ( params.x(idx+1) - params.x(idx) ), 0];

		if isfield(params,'membrane_func')
			% Set Membrane Capacitance
			Nm = params.membrane_N;
			Cspecific = 1e-2; % Membrane capacitance in farads per m^2
			params.Cl(Nm) = Cspecific * params.A(Nm+1) ;
		end

           end

        % Transverse capacitance
	   if ~isfield(params,'Ct');

		if ~isfield(params,'d_ratio')
			params.d_ratio = 0.75; % Inner pipette diameter / Outer Pipette Diameter
		end

		% Position of boundaries between each section
		b = ( params.x(idx) .* R(idx+1) + params.x(idx+1) .* R(idx) ) ./  ( R(idx) + R(idx+1) );	

		% Cone angle between each section.
		theta = 2 * atan(  ( R(idx+1)- R(idx) ) ./ (params.x(idx+1) - params.x(idx))  );
		if isfield(params,'membrane');	theta(params.membrane_N) = 0; end

		% Transverse capacitance of each element.
		epsilon_wall = 3 * 8.85e-12; % permittivity of walls
		Ctb = [0, 2 * pi * epsilon_wall  ./ log(1./params.d_ratio) .* (params.x(idx+1) - b) ];
		Ctf = [2 * pi * epsilon_wall ./ log(1/params.d_ratio) .* (b - params.x(idx)), 0];
                
		% Whole Cell - Adjust for bath
		if strcmp(params.config,'whole-cell')
			Ctf(params.Nseg-1) = 0;
			Ctb(params.Nseg) = 0;
		end
		params.Ct = Ctb + Ctf;
       end	
       
	% Solution parameters
	% If not specified, assume 150mM KCl
	  % Diffusion Coefficients
	   if ~isfield(params,'D'); 
		params.D = 1.957e-9* [1 1.038]';      % Diffusion coefficient of K+,  Cl-  in m^2/s
						   % http://web.med.unsw.edu.au/phbsoft/mobility_listings.htm
					           % http://pubs.acs.org/doi/pdf/10.1021/ja01622a016
	   end	

          % Fundamental Charge per molecule
	    if ~isfield(params,'z')
		params.z = [1 -1]'; % Default is KCl
	    end

  	  % Solution in Pipette
       if ~isfield(params,'rho_pip')
    		params.rho_pip = [150 150]'; % Concentration in pipette - 150mM KCl
	   end    
        
	  % Solution in Bath    
	   if ~isfield(params, 'rho_bath')
	 	    params.rho_bath = params.rho_pip;
       end
      
	% Concentration at each reference point
	if ~isfield(params,'rho')
            if isfield(params,'membrane_func')
			% Membrane Present
			% Concentration params.rho_pip on pipette side; concentration params.rho_bath on bath side
        	    	if strcmp(params.config,'inside-out');
		                sf = ([1:params.Nseg] > params.membrane_N);
		            else
        		        sf = ([1:params.Nseg] <= params.membrane_N);
	        	    end
	    	else
			% No membrane present
			% Initialize as linear gradient
                    sf = 1 - ([1:params.Nseg]-1) / (params.Nseg - 1);
            end
		
		    % Determine bath concentrations
           	if isa(params.rho_bath,'function_handle')
               		rho_b = params.rho_bath(params.t);
           	else
               		rho_b = params.rho_bath;
           	end
        
	        % Set concentrations
            params.rho = params.rho_pip * sf + rho_b * (1-sf) ;
        end

	% Holding voltage
	    if ~isfield(params,'V_hold')
        	params.V_hold = 0;
    	    end
   
    	% Location of pipette and bath electrodes
    	   if (~isfield(params,'n_pip'))|(~isfield(params,'n_bath'))
          	if strcmp(params.config,'inside-out')
               		params.n_pip = params.Nseg;
                	params.n_bath = 1;
            	else
                	params.n_pip = 1;
                	params.n_bath = params.Nseg;
            	end
           end
              
      % Voltage at each point in system
    	   if ~isfield(params,'V')
         	params.V = zeros(1,params.Nseg);
	        params.V(params.n_pip) = params.V_hold * (1 - 2 *(params.n_pip > 1));
        	params.V = voltage_gauss(params); 
           end
          
      % Fluxes
	   if ~isfield(params,'J')  % Molecular fluxes 
	 	params.J = zeros(length(params.z), params.Nseg-1);
	   end

           if ~isfield(params,'Jsol') % solvent flux
    		params.Jsol = zeros(1,  params.Nseg-1);
           end
    
           % Creation/Conversion Rates
	   if ~isfield(params,'rr')
		params.rr = zeros(length(params.z), params.Nseg); % Reaction Rates
       end

       % Series Resistance compensation
       % If Rseries is not defined, no series resistance compensation.
       if isfield(params, 'Rseries')
           if ~isfield(params,'Rseries_tau');
               params.Rseries_tau = 1e-5; % Time constant in seconds
                                          % Stable when tau > Cpipette * Rpip                                        % R
           end
       end
       
       % Define constraint functions.
       if ~isfield(params, 'Constraints')
               params.Constraints.type = 'basic';   % Basic constraints on just rho, J and V
       end		
            % Constraints.type -> 	optional name describing constraints  (e.g. 'basic, 'buffered', 'solvent', etc.)
            % Constraints.param     -> list of parameter names being constrained.
            % Constraints.fhandle   -> function constraining each of listed parameters
            % Constraints.segments  -> compartments to which constraints are applied
            % Constraints.components-> components of parameter to which constraints are applied
            % Constraints.stepsize  -> amount to vary parameter by when calculating derivatives
            % Constraints.err_tol    -> residual error allowed after optimization.
       if ~isfield(params.Constraints,'param');

           switch  params.Constraints.type

               	case {'basic'}
				% Adjust concentrations (rho), molecular fluxes (J), and voltage (V)
				params.Constraints.param =      {'rho',             'J',                    'V'}; 
				params.Constraints.fhandle =  	{@constraint_Rho,	@constraint_J,          @constraint_V};
                params.Constraints.segments = { [2,params.Nseg-1],      [1,params.Nseg-1],      [1,params.Nseg],};
    			params.Constraints.components ={[1,size(params.rho,1)], [1,size(params.rho,1)], [1 1]};
    			params.Constraints.stepsize = 	{0.001,                   0.01,                 1e-5};
   				%params.Constraints.err_tol = 	{0.001,                     1e-18,              1e-5};
                params.Constraints.err_tol = 	{0.0000001,                     1e-18,              1e-5};

                               
                case {'solvent'}
				% Adjust concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol)
				params.Constraints.param =	{'rho', 		'J',                    'V',			'Jsol'}; 
				params.Constraints.fhandle =  	{@constraint_Rho,	@constraint_J,          @constraint_V, 		@constraint_Jsol};
    			params.Constraints.segments = { [2,params.Nseg-1],      [1,params.Nseg-1],      [1,params.Nseg],	[1,params.Nseg-1]};
    			params.Constraints.components ={[1,size(params.rho,1)], [1,size(params.rho,1)], [1 1],			[1 1]};
    			params.Constraints.stepsize = 	{0.001,                   0.01,                 1e-5,			1e-18};
   				params.Constraints.err_tol = 	{0.001,                     1e-18,              1e-5,			1e-19};

                case {'buffer'}
				% Adjust concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
				params.Constraints.param =	{'rho', 		'J',                    'V',			'Jsol',			'rr'}; 
				params.Constraints.fhandle =  	{@constraint_Rho,	@constraint_J,          @constraint_V, 		@constraint_Jsol, 	@constraint_Chemistry_AB};
    			params.Constraints.segments = { [2,params.Nseg-1],      [1,params.Nseg-1],      [1,params.Nseg],	[1,params.Nseg-1],	[2, params.Nseg-1]};
    			params.Constraints.components ={[1,size(params.rho,1)], [1,size(params.rho,1)], [1 1],			[1 1],			[1 3]};
    			params.Constraints.stepsize = 	{0.001,                   0.01,                 1e-5,			1e-18,			1};
                params.Constraints.err_tol = 	{0.001,                     1e-18,              1e-5,			1e-19,			1e-3};

                case {'basic_perfusion'}
				% Adjust concentrations (rho), molecular fluxes (J), and voltage (V)
                % Allow concentrations of pipette and bath compartments to
                % evolve in response to perfusion.
				params.Constraints.param =      {'rho',             'J',                    'V'}; 
				params.Constraints.fhandle =  	{@constraint_Rho,	@constraint_J,          @constraint_V};
                params.Constraints.segments = { [1,params.Nseg],      [1,params.Nseg-1],      [1,params.Nseg],};
    			params.Constraints.components ={[1,size(params.rho,1)], [1,size(params.rho,1)], [1 1]};
    			params.Constraints.stepsize = 	{0.001,                   0.01,                 1e-5};
   				params.Constraints.err_tol = 	{0.001,                     1e-18,              1e-5};
                
                case {'solvent_cell_volume'}
				% Water flow across whole-cell membrane, but no water flow in pipette
                % Adjust concentrations (rho), molecular fluxes (J), voltage (V), solvent flux at membrane (Jsol), and cell volume
                idx_cell = [params.membrane_N,params.membrane_N]; % Only consider water flux at membrane and changes in cell volume
				params.Constraints.param =      {'rho',                 'J',                    'V',			'Jsol',             'Lambda'}; 
				params.Constraints.fhandle =  	{@constraint_Rho,	@constraint_J,          @constraint_V, 		@constraint_Jsol,   @constraint_Lambda};
    			params.Constraints.segments = { [2,params.Nseg-1],      [1,params.Nseg-1],      [1,params.Nseg], idx_cell,          idx_cell };
    			params.Constraints.components ={[1,size(params.rho,1)], [1,size(params.rho,1)], [1 1],			[1 1],              [1,1]};
    			params.Constraints.stepsize = 	{0.001,                   0.01,                 1e-5,			1e-18,              1e-20};
   				params.Constraints.err_tol = 	{0.001,                     1e-18,              1e-5,			1e-19,              1e-20};
                
                
		    end		
       end
      
	% Calculate optimal ordering of the constrained variables
	% Solutions are determined by finding values of the N_constraints constrained variables 
	% (each described by the parameter j, compartment n, component m) so each satisfies its constraint equation.
	%
	% Because the constraints are all local (i.e. the value of any constraint function 
	% depends only on the value of variables in adjacent compartments), there are several advantages
	% to ordering the constrained variables by compartment.  
	%	1. The Hessian matrix is banded
	%	2. The Hessian can be evaluated using far fewer function calls.
	%	
	% However, the constraint functions are most naturally ordered and evaluated by parameter,
	% 
	% Thus, we want to define indices to allow easy transposition of variables between ordering by parameter and ordering by compartment.
	% To do this, for each parameter we can define a pair of indices 
	% 	idx_v =  params.Constraints.idx_p{j} -> list of all variables in parameter which are constrained.
	%						e.g. [2, 3, ..., Nseg-1]
	% 	idx_c =  params.Constraints.idx_c{j} -> ordering of these variables by compartment.
	% 
	% By making the idx_v(k)-th variable of parameter j to be the idx_c(k)-th constraint, the constrained variables are then
	% listed by compartment.
	% Results can then be easily mapped in either direction.
	% To evaluate constraints, the constraint function, cfunc_j, of parameter j then maps
	% to constraint errors, 
	% 	i.e. error(idx_c) = cfunc_j(idx_v)
	% Conversely, new values for variables ordered by constraint, vari, can be mapped back to parameters by
	%	     params_j(idx_v) = vari(idx_c)

    	  % Calculate number of constraints
    	 	  for j=1:length(params.Constraints.param) % Step through each parameter
        		sidx = params.Constraints.segments{j}; 	% range of compartments which are constrained
                if (length(sidx)==1);   sidx(2) =sidx(1);  end; 
                cidx = params.Constraints.components{j};% range of components which are constrained
        		Nc(j) = (sidx(2)-sidx(1)+1) * (cidx(2)-cidx(1)+1); % Number of variables constrained for this parameter
	   	  end
	   	  N_constraints = sum( Nc); % Number of constraints
		  params.Constraints.N = N_constraints;
    
    	  % Now list all constraints in terms of output order, parameter, compartment and component
    	   	  idx_ord = [1:N_constraints];       % Output order
    	   	  idx_para = zeros(1,N_constraints); % Index of variable being constrained (within the parameter)
    	   	  idx_seg = zeros(1,N_constraints);  % Constrained variable is in compartment idx_seg and component idx_comp
    	   	  idx_comp = zeros(1,N_constraints); 

 	          ctr = 0;
   		  for j=1:length(params.Constraints.param) 	 % Step through each constrained parameter in turn

        		sidx = params.Constraints.segments{j};	 % Range of compartments which are constrained 
                if (length(sidx)==1);   sidx(2) =sidx(1);  end;
        		cidx = params.Constraints.components{j}; % Range of components which are constrained 
        		idx_loc = ctr + [1:Nc(j)];		 % Index of constraints
		
                % For each element, define parameter, compartment and component being constrained.
        		idx_para(idx_loc) = j;			 
                idx_seg(idx_loc)  = reshape( repmat([sidx(1):sidx(2)],[(cidx(2)-cidx(1)+1),1]),[1,Nc(j)]); 
                idx_comp(idx_loc) = reshape( repmat([cidx(1):cidx(2)],[1,(sidx(2)-sidx(1)+1)]),[1,Nc(j)]); 
                    
        		ctr = ctr + Nc(j);
    	   	  end   

    	  % Sort constraints by compartment number
   	       		[idx_seg, idx_r] = sort(idx_seg);
    			idx_ord = idx_ord(idx_r);   % Index of output ordered by compartment
    			idx_para = idx_para(idx_r);   % Index of parameters ordered by compartment
    			idx_comp = idx_comp(idx_r); % Index of components ordered by compartment
                
	% Organize this mapping by parameter  
      params.Constraints.idx_c = cell(1,length(params.Constraints.param));
      params.Constraints.idx_v = cell(1,length(params.Constraints.param));
  	  for j=1:length(params.Constraints.param)
            sidx = params.Constraints.segments{j};	
            idx_t = find(idx_para==j);
        	% Index of each constrained variable in constraint space
            params.Constraints.idx_c{j} = idx_t; 
            
            % Index of each constraint variable in parameter space
            if length(sidx)==2 
                % Constraint function applies to multiple compartments 
                params.Constraints.idx_v{j} = idx_comp(idx_t) + (idx_seg(idx_t)-1) * size(params.(params.Constraints.param{j}),1); 
            else
                % Constraint has only a single compartment and only has a
                % component index
                params.Constraints.idx_v{j} = idx_comp(idx_t) ; 
            end
      end

function result = voltage_gauss(params)
% 	Apply Gauss's law to determine voltage at every point in system
% 	Assume voltages are fixed in first and last compartments
%
% 	Input
%		params -> state of system
%
% 	Output 
%		result = voltage at each point

	Np = params.Nseg;

    if (Np >2)
        % Determine limits
        Na = 2; Nb = Np-1;
        idx = [Na:Nb];
    
        % Charge
        Q = sum((params.F * params.z * params.Lambda(idx) ) .* params.rho(:,idx) , 1);

        % Capacitances at clamped ends.
        if (Na>1) ;  Q(1) = Q(1) + params.Cl(1) * params.V(1); end; % Charge from left side
        if (Nb<Np) ; Q(1+Nb-Na) = Q(1+Nb-Na) + params.Cl(Np-1) * params.V(Np); end; % Charge from right side

        % Capacitance
        Clt = [0, params.Cl, 0];
        Cl = params.Cl([Na:(Nb-1)]);
        C = diag( params.Ct(idx)+ Clt([Na:Nb]) + Clt(1+[Na:Nb]) ) - diag(Cl,1) - diag(Cl,-1); 
     
    	% Solve for voltage
    	result = params.V;
    	result(idx) = Q/C;
    else
        result = params.V;
    end

    