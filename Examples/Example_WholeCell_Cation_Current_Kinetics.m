% function Example_WholeCell_Cation_Current
%        Calculate evolution of I and rho versus V for step change in voltage.
%           150mM KCl in pipette and bath, 
%           Raccess = 5 MOhms.
%           g_K+ = 20 nS

function Example_WholeCell_Cation_Current_Kinetics

% Configuration parameters
		params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                   
	% Solution Parameters
        params.D =  [1,  1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ and Cl-.
		params.z =  [1, -1]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [150,  150  ]';  	% 150mM NaCl in bath
		params.rho_pip =  [150,  150  ]';  	% 150mM KCl in pipette
        
	% Pipette parameters
        R_pip = 5e6 ; % Pipette resistance of 5 MOhms
		params.theta = 10 * pi/180;  % Tip Angle
        params.r_tip = calculate_tip_diameter(params, R_pip)/2;% Tip radius of pipette (in m)
		                            
	% Cell and Membrane parameters
        params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0]'; 		% Initial Membrane permability to each ion in m/s
        params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
		params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell - 1 pL = 1e-15 m^3
        
	% Initialize all remaining parameters
        params = initialize_system(params);
    
    % Calculate membrane permeability for an initial current of 2nA K+ at % 100mV (i.e. 
        params.V_hold = 100e-3; 	% Initial holding voltage
        I_mem = [2e-9, 0 ]'; % Initial current of 2nA by K+
        P_GHK = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, params.V_hold, params.z, params.A_cell, params.T);

    % Calculate Time-Dependent Evolution.
        params = evolve_backward_euler(params,inf); % Initialize parameters with g_mem = 0nS
        params.membrane_perm = P_GHK; % Set membrane permeability.
        t_values = [0.001, 0.005,0.1, 0.5, 1:1:20]; % Time points to evalute from 0 to 20 second interval
        
        [params, results] = evolve_backward_euler(params, t_values); % Calculate SS 
      
      % Plot Current and concentration versus t.
        figure(1); clf;
            subplot(2,1,1); % Current versus time
                plot(t_values, results.I_pip * 1e9,'k');
                ylabel('I (nA)');
                title('Step Increase in Membrane K+ Conductance');
            subplot(2,1,2); % Concentration versus time
                plot(t_values, results.rho_intra(1,:),'k');
                ylabel('\rho_{cell} (mM)');
                xlabel('t (s)');
        
     