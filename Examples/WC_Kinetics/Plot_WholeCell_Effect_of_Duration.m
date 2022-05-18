% function plot_WholeCell_Effect_of_Duration
%
%   Plot Current flow during voltage ladder
%
%   


function plot_WholeCell_KCl_ladder

	% Configuration parameters
		params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
	% Solution Parameters
		
        params.D =  [1,  1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ and Cl-.
		params.z =  [1, -1]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [150,  150  ]';  	% 150mM NaCl in bath
		params.rho_pip =  [150,  150  ]';  	% 150mM KCl in pipette
        
	% Pipette parameters
        R_pip = 5e6 ; % Pipette resistance of 5 MOhms
		params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
        params.theta = 10 * pi/180;  % Tip Angle
        [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip); % Calculate tip radius
        params.r_tip = tip_diameter/2;% Tip radius of pipette (in m)
		params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
        params.Nseg  =  11; 		% Number of segments between the pipette electrode and tip
		params.V_hold = 0e-3; 	% Initial holding voltage
        params.Rseries = 4.5e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
                                    
	% Cell and Membrane parameters
        params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0]'; 		% Membrane permability to each ion in m/s
        params.membrane_Pf = 0e-5   ; 			% Membrane permeability to water (m/s)
		params.A_cell = 1e-9;                  % Cell Area needed for whole-cell
		params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell

	% Initialize all remaining parameters
    params = initialize_system(params);
    
    % Calculate membrane permeability for an initial current of 2nA K+
    I_mem = [10e-9, 0 ]'; % Initial current of 10nA by K+
    V_mem = 100e-3; % 10nA for 100mV 
    P_GHK = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, V_mem, params.z, params.A_cell, params.T);
    
    % Report Access Resistance and Series Resistance Correction
   	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(params)/1e6), params.Rseries/1e6)
    fprintf('\n Membrane Permeability = [%e, %e] m/s', P_GHK);
    
    
    paramst = params;
    
    V = [-100:25:100] * 1e-3;
    
    dval = [50e-3, 1e-5, 9.9e-4, 9e-3,   9e-2,   5e-3,    45e-3 ] ; % Duration of each of intervals
    intervals.N_eval =   [2,      10,   10,     10,     10,     10,     20] ; % Number of timepoints
    
    figure(1); clf;

    
    params.membrane_perm = P_GHK;
    params.Ct = 0 * paramst.Cl;
    params.Cl = 0 * paramst.Cl;
    paramst = params;
    paramst.V_hold = 0;
    [paramst, results] = evolve_backward_euler(paramst, inf);
    
    sf = [1, 20]; % Scale time
    
    for k=1:2
    
        for j=1:length(V)

        intervals.variables.V_hold = {0e-3, V(j), V(j), V(j), V(j), 0e-3, 0e-3};
        intervals.duration =  dval * sf(k);
        paramst = params;
        protocol = create_protocol( intervals );
        [paramst, results] = evolve_backward_euler(paramst, protocol);
        
        paramst.V_hold = V(j);
        [paramst, results_f] = evolve_backward_euler(paramst, inf);
        
        subplot(2,2,k);
        plot(protocol.t, results.I_pip*1e9,'k'); hold on;

        subplot(2,2,k+2);
            plot(protocol.t, mean(results.rho_intra),'k'); hold on;
        
        end
    end
    
        
        subplot(2,2,1)
            axis([0 0.2 -11.2 11.2])
            xlabel('t (s)');
            ylabel('I (nA)');
            plot( [0.01 0.01], [0 5],'k','Linewidth',1); % 5nA scale bar
            plot( [0.0 0.05], [1, 1],'k','Linewidth',1); % 50ms scale bar
            axis off
            
        subplot(2,2,3)
            xlabel('t (s)')
            ylabel('\rho (mM)');
            axis([0 0.2 144 156])
            box off
            set(gca,'YTick',[145 150 155],'TickDir','out');
            
            
        subplot(2,2,2)
            axis([0 4 -11.2 11.2])
            xlabel('t (s)');
            ylabel('I (nA)');
            plot( [0.1 0.1], [0 5],'k','Linewidth',1); % 5nA scale bar
            plot( [0 1], [1, 1],'k','Linewidth',1); % 1s scale bar
            axis off
            
            
        subplot(2,2,4)
            xlabel('t (s)')
            ylabel('\rho (mM)');
            axis([0 4 55 245])
            box off;
            set(gca,'YTick',[100 125 150 175 200 225],'TickDir','out');
    
	save_figure('fig_WholeCell_Effect_of_Duration',[ 15 15]);       
            
    