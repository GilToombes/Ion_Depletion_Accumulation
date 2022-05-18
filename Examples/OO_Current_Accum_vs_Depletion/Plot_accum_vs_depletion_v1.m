% function fig_accum_vs_depletions_v1
%
% Illustrate different impact of accumulation and depletion on current
% Outside-out membrane patch
%     R = 5 MOhms
%     150mM KCl
%

function fig_accum_vs_depletions_v1

	% Configuration parameters
		params.config = 'outside-out';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
	% Solution Parameters
        params.D =  [1,  1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ and Cl-.
		params.z =  [1, -1]';	  	% Relative charge of K+ and Cl-.
		params.rho_bath = [150, 150  ]';  	% 150mM KCl in bath
		params.rho_pip =  [150, 150  ]';  	% 150mM NaCl in pipette
        
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
        params.membrane_x = 1e-9;   			% Position of membrane patch back from the tip (in m)
		params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0]'; 		% Membrane permability to each ion in m/s

	% Initialize all remaining parameters
    params = initialize_system(params);
    
 % Calculate membrane permeability for an initial current of 2nA K+
    I_mem = [2e-9, 0 ]'; % Initial current of 2nA by K+
    V_mem = 100e-3; % Voltage in mV
    P_GHK = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, V_mem, params.z, params.A(params.membrane_N), params.T);
    
    % Report Access Resistance and Series Resistance Correction
   	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(params)/1e6), params.Rseries/1e6)
    fprintf('\n Membrane Permeability = [%e, %e] m/s', P_GHK); % Calculate membrane permeability for an initial current of 2nA K+
    
    paramst = params;
    
    V = [-100:25:-25,  0.001,  25:25:100] * 1e-3;
    
    intervals.duration =  [25e-3, 1e-5, 9.9e-4, 9e-3, 9e-2  5e-3, 20e-3 ] ; % Duration of each of intervals
    intervals.N_eval =   [2,    40, 10,     10, 10,     10, 20] ; % Number of timepoints
    
    figure(1); clf;
                    % Left,     Bottom,     Width,  Height
        panel_pos = [ 0.05,     0.7,        0.5,    0.25 ;
                      0.05,     0.4,        0.5,    0.25 ;
                      0.65,     0.7,        0.3,    0.25 ;
                      0.65,     0.4,        0.3,    0.25;
                      0.05,     0.05,        0.5,   0.15;
                      0.65,     0.05,         0.3,   0.25; ];
                      
                      
    params.membrane_perm = P_GHK;
    paramst = params;
    paramst.V_hold = 0;
    [paramst, results] = evolve_backward_euler(paramst, inf);
    
    for j=1:length(V)
        
        %paramst = params;
        %paramst.V_hold = V(j);
        %[paramst, results] = evolve_backward_euler(paramst, inf);
        %paramst.membrane_perm = 1.1e-4 * [1 0]';
        
        intervals.variables.V_hold = {0e-3, V(j), V(j), V(j), V(j) ,0e-3, 0e-3};
        paramst = params;
        protocol = create_protocol( intervals );
        [paramst, results] = evolve_backward_euler(paramst, protocol);
        
        paramst.V_hold = V(j);
        [paramst, results_f] = evolve_backward_euler(paramst, inf);
        
        subplot('Position',panel_pos(1,:));
        plot(protocol.t, results.I_pip*1e9,'k'); hold on;

        subplot('Position',panel_pos(2,:));
            plot(protocol.t, mean(results.rho_intra),'k'); hold on;
        
        I_max(j) = max(abs(results.I_pip)) * sign(results_f.I_pip);
        I_fin(j) = results_f.I_pip;
        rho_fin(j) = mean(results_f.rho_intra);
        
    end
    
        subplot('Position',panel_pos(1,:))
            axis([0 0.15 -2.1 2.1])
            plot([0.002, 0.002], [0.1, 1.1],'k'); % Add 1nA scale bar
            axis off;
       
        subplot('Position',panel_pos(2,:))
            axis([0 0.15 95 205])
            plot([0.002 0.052], [100, 100], 'k'); % Add a 50ms time bar
            plot([0.002 0.002], [160, 200], 'k'); % Add a 40mM scale bar
            axis off
            
        subplot('Position',panel_pos(3,:))
            plot(V*1e3, I_max*1e9,'k.-'); hold on
            plot(V*1e3, I_GHK(paramst, V, 0,150)*1e9,'k--'); hold on;
            plot(V*1e3, I_GHK(paramst, V, 150, 0)*1e9,'k--'); 
            
            xlabel('V (mV)');
            ylabel('I (nA)');
            axis([-120 120 -2.2 2.2])
            ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
            set(gca,'XTick',[-100 -50 50 100]);
            set(gca,'YTick',[-2 -1 1 2]);
           
        subplot('Position',panel_pos(4,:))
            plot(V*1e3, I_fin*1e9,'k.-'); hold on
            plot(V*1e3, I_GHK(paramst, V, rho_fin, 0)*1e9, 'k--');
            plot(V*1e3, I_GHK(paramst, V, 0,150)*1e9,'k--');
            
            axis([-120 120 -2.2 2.2])
            xlabel('V (mV)');
            ylabel('I (nA)');
            set(gca,'XTick',[-100 -50 50 100]);
            set(gca,'YTick',[-2 -1 1 2]);
            ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
   
        % Now Measure Vrev at end of puls
        paramst = params;
        V = [-15:5:15] * 1e-3;
    
        intervals.duration =  [25e-3,   1e-5,   9e-5,     9e-4,   9e-3, 9e-2  1e-5,   9e-5,   9e-4,   9e-3,   9e-5, 9e-4, 14e-3 ] ; % Duration of each of intervals
        intervals.N_eval =    [2,       10,     10,         10,     10, 10,    10,     10,     10,       10,   10, 10,15    ] ; % Number of timepoints
        idx_ladder =    2 + 10*5 + 10* 2 ; % index for 1msec after step into voltage ladder.
        
        paramst = params;
        paramst.V_hold = 0;
        [paramst, results] = evolve_backward_euler(paramst, inf);

        for j=1:length(V)
        
        %paramst = params;
        %paramst.V_hold = V(j);
        %[paramst, results] = evolve_backward_euler(paramst, inf);
        %paramst.membrane_perm = 1.1e-4 * [1 0]';
        Vstep = -100e-3;
        intervals.variables.V_hold = {0e-3, Vstep, Vstep,   Vstep,  Vstep,  Vstep, V(j),   V(j),   V(j),   V(j),   0e-3, 0e-3, 0e-3};
        paramst = params;
        protocol = create_protocol( intervals );
        [paramst, results] = evolve_backward_euler(paramst, protocol);
        
        paramst.V_hold = V(j);
        [paramst, results_f] = evolve_backward_euler(paramst, inf);
        
        subplot('Position',panel_pos(5,:));
        plot(protocol.t, results.I_pip*1e9,'k'); hold on;
       
        %subplot('Position',panel_pos(6,:));
        %plot(protocol.t, protocol.variables.V_hold *1e3,'k'); hold on;
        
        I_ladder(j) =    results.I_pip(idx_ladder);
        V_ladder(j) =  V(j);
        
        end
        
       
        subplot('Position',panel_pos(5,:))
            axis([0 0.15 -2.1 0.5])
            plot([0.002, 0.002], -1 * [0.1, 1.1],'k'); % Add 1nA scale bar
            axis off;
       
        subplot('Position',panel_pos(6,:));
            plot(V_ladder*1e3, I_ladder*1e9, 'k.-');
            ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
            axis([-15 15 -0.3 0.6]);
            set(gca,'XTick',[-10 -5 5 10]);
            xlabel('V (mV)');
            ylabel('I (nA)');
            
    	save_figure('Fig_accum_vs_depletion_v1',[ 15 20]);       

    
    function result = I_GHK(params, V, rho_in, rho_out)
        
        idx = 1;
        Am = params.A(params.membrane_N);
        u = params.z(idx)*params.F*V/(params.R*params.T) ; 
        uz = (abs(u)<1e-6); % Uncharged particles and zero membrane voltge
        result = params.F * Am * params.membrane_perm(idx) .* (u + uz) ./ (1 - exp(-u) + uz) .* (rho_in - rho_out.*exp(-u));