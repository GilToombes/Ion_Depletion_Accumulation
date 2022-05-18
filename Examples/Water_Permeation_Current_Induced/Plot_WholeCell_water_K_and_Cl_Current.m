% function Fig_WholeCell_Water_K_Current


function Fig_WholeCell_Water_K_and_Currents

	% Configuration parameters
        params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'solvent'; % Solvent flow through both membrane and pipette
    	
	% Solution Parameters
        params.D =  [1, 0.33, 0.682, 1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ , gluconate and Cl-.
		params.z =  [1, -1, 1 ,  -1]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [0,     0, 150, 150 ]';  	% 150mM NaCl in bath
		params.rho_pip =  [150,150 , 0,   0   ]';  	% 150mM KCl in pipette
        
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
		params.membrane_perm =  [ 0 0 0 0]' ; 		% Membrane permability to each ion in m/s
        params.membrane_Pf = 5e-5  ; 			% Membrane permeability to water (m/s)
	    params.A_cell = 1e-9;                  % Cell Area needed for whole-cell   Boss et. al. 2013
	    params.Lambda_Cell = 1e-15 ;           % Cell volume needed for whole-cell
        
        params = initialize_system(params);
        [params, results] = evolve_backward_euler(params, inf);
        
        
        % Membrane Permeabilities
        
            % K Current
            I_mem = [2e-9, 0, 0,0 ]'; % Initial current of 1nA by K+
            V_mem = params.V_hold ; % Voltage in mV
            P_GHK_K = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, V_mem, params.z, params.A_cell, params.T);
        
            % Cl Current
            I_mem = [0, 0, 0,2e-9 ]'; % Initial current of 1nA by K+
            V_mem = params.V_hold ; % Voltage in mV
            P_GHK_Cl = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, V_mem, params.z, params.A_cell, params.T);
        
            Pf = 5e-5; % Membrane Water Permeability
            
        % Report Access Resistance and Series Resistance Correction
        fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
   	    fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(params)/1e6), sum(params.Rseries/1e6))
        fprintf('\n Membrane Ionic Permeability = [%e, %e, %e, %e] m/s', P_GHK_K); % Membrane ionic Permeability
        fprintf('\n Membrane Ionic Permeability = [%e, %e, %e, %e] m/s', P_GHK_Cl); % Membrane ionic Permeability
        fprintf('\n Membrane Water Permeability = %e m/s', Pf); 
        
        intervals.duration=                 [2,         0.1,        1.9,        28,     0.1,        1.9     ]; % Duration of each of intervals
        intervals.N_eval =                  [2,         5,          10,          26,     5,        20       ]; % Number of timepoints during each interval
        
        P_GHK = P_GHK_K; 
        intervals.variables.membrane_perm=  {0*P_GHK,   P_GHK,    P_GHK,   P_GHK,    0*P_GHK,    0*P_GHK };% Membrane Permeability.
        protocol_K = create_protocol(intervals); % Create Protocol 
        
        P_GHK = P_GHK_Cl; 
        intervals.variables.membrane_perm=  {0*P_GHK,   P_GHK,    P_GHK,   P_GHK,    0*P_GHK,    0*P_GHK };% Membrane Permeability.
        protocol_Cl = create_protocol(intervals); % Create Protocol 

        figure(1); clf;
        
            % Panel Positions
            xpos =[ 0.1, 0.6];
            width =  [0.35, 0.35];
            ypos = [ 0.72, 0.52, 0.26, 0.05];
            height = [0.25, 0.18, 0.18, 0.15];
            
            protocol = {protocol_K, protocol_Cl};
            idx_ion = [1, 4]; % K+ in first protocol, Cl- in second
            
            for k=1:length(protocol)
                
                paramst = params; % No water permeation
                paramst.membrane_Pf = 0e-5  ; 
                [paramst, results_nowater, params_full_nowater] = evolve_backward_euler(paramst, protocol{k});
                
                paramst = params;   % Water permeation
                paramst.membrane_Pf = 5e-5  ; 
                [paramst, results_water, params_full_water] = evolve_backward_euler(paramst, protocol{k});
                
                
                % Plot results
                t_val = protocol{k}.t;
                t_range = [0 34];
                
                % Current
                subplot('Position',[xpos(k), ypos(1), width(k), height(1)] ); 
                        plot(t_val, results_water.I_pip*1e9,'k'); hold on
                        plot(t_val, results_nowater.I_pip * 1e9, 'k--');
                 
                        
                        plot( 6 + [0 10], 0.3 + [0 0 ], 'k'); % 10 second time bar
                        plot( 1 + [0 0], 0.3 + [0 1 ], 'k'); % 1nA current bar
                        axis([ t_range -0.1 2.1]);
                        axis off
                
                % Concentration          
                subplot('Position',[xpos(k), ypos(2), width(k), height(2)] ); 
                        plot(t_val,   results_water.rho_intra(idx_ion(k),:),'k'); hold on    
                        plot(t_val, results_nowater.rho_intra(idx_ion(k),:), 'k--');
                        switch k
                            case 1
                                axis([ t_range  100 152]);                              
                                set(gca,'YTick',[110 130 150]);
                                ylabel('\rho_{cell,K+} (mM)');
                            case 2
                                axis([ t_range  0 52]);         
                                set(gca,'YTick',[0 20 40]);
                                ylabel('\rho_{cell,Cl-} (mM)');
                        end
                        box off;
                        h=gca; h.XAxis.Visible='off';
                        
            % Flow
            subplot('Position',[xpos(k), ypos(3), width(k), height(3)] ); 
                for j=1:length(params_full_water)
                    Jsol(j) =           params_full_water(j).Jsol(params.membrane_N);
                    Jsol_nw(j) =           params_full_nowater(j).Jsol(params.membrane_N);
                    %Lambda_cell(j) =    params_full(j).Lambda(params{k}.membrane_N);
                    %IK(j) =             params_full(j).J(1,params{k}.membrane_N-1) * params{k}.F;
                    %ICl(j) =            params_full(j).J(2,params{k}.membrane_N-1) * -1 * params{k}.F;
                end
            
                    plot(t_val, Jsol*1e18,'k'); hold on
                    plot(t_val, Jsol_nw*1e18,'k:'); hold on
                    
                    axis([ t_range  -30 30]);                              
                    box off;
                    h=gca; h.XAxis.Visible='off';
                    ylabel('f_{mem} (fL.s^{-1})');
                    set(gca,'YTick',[-20 0 20]);    
           
                 
                        
            % Volume
            %subplot('Position',[xpos(k), ypos(4), width(k), height(4)] ); 
            %        plot(t_val, Lambda_cell*1e15,'k'); hold on
                    %    axis([  t_range 0.6 1.02]);
             %           box off;
             %           h=gca; h.XAxis.Visible='off';
             %           ylabel('\Lambda_{cell} (pL)');
             %           set(gca,'YTick',[0.6 0.8 1]);    
            end
        
    save_figure('Fig_Whole_Cell_Water_K_and_Cl_Current_v1',[8 6.5])