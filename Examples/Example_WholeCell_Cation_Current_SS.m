% function Example_WholeCell_Cation_Current_SS
%
%       Calculate Steady-state I and rho versus V
%           150mM KCl in pipette and bath, 
%           Raccess = 5 MOhms.
%           g_K+ = 20 nS

function Example_WholeCell_Cation_Current_SS

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
        params.membrane_perm = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, params.V_hold, params.z, params.A_cell, params.T);
    
    % Calculate Steady-state I and rho versus V
        V_hold = [-100, -75, -50, -25, 0, 25, 50, 75, 100] * 1e-3; % Holding voltages in mV
        
        for j=1:length(V_hold) % Calulate for each voltage in turn.
        
             params.V_hold = V_hold(j);    % Set Holding Voltage
            [params, results] = evolve_backward_euler(params, inf); % Calculate SS 
            I_SS(j) = results.I_pip ;    % SS Current
            rho_SS(j) = results.rho_intra(1); % SS K+ Concentration 
        end
    
    % Plot I versus V
        figure(1); clf;    
            plot(V_hold * 1e3 , V_hold * 20e-9 * 1e9, 'r'); hold on; % Current for ideal 20nS conductance
            plot(V_hold * 1e3, I_SS * 1e9, 'k'); % Steady-State Current
            xlabel('V_{hold} (mV)');
            ylabel('I (nA)');
            legend('Ideal','Steady-State','Location','SouthEast');
            ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; box off % Cross axes at origin

    % Plot rho versus I
        figure(2); clf;
            plot(I_SS * 1e9, rho_SS , 'k');
            xlabel('I (nA)');
            ylabel('\rho_{cell} (mM)');
    
     