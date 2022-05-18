% function Figure_9_Tau_Pipette_v2
%
%   Show effect of tip diameter and angle on time dependence of ion
%   depletion
%   Example has 5 MOhm pipette containing 150mM KCl
%   Outside-out config with K-channels giving an initial current of 1nA (i.e. 5mV voltage-drop)
%   Calculate time-course for 3 pipettes, each with 5 MOhm ressitance but
%           theta_tip - 6.66, 10, and 15 degrees
%           d_tip       0.96, 0.64, and 0.43 microns
%   Show concentration and then concentration scaled to
%                       tau_diff = (d_tip / (2*sin(theta/2)))^2 / D
%   
%   Then show plots of  R and tau_diff versus d_tip and theta_tip

function Figure_Time_Dependence_dtip_theta

	% Configuration parameters
		params.config = 'outside-out';      % Options - 'inside-out','outside-out','whole-cell'
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
                                            %       solvent - voltage, molecular concentrations and currents + solvent flow
                                            %       buffer  - concentrations (rho), molecular fluxes (J), voltage (V), and solvent flux (Jsol), and buffered rr
                                            % Alternatively, specify a set of constraints.
	% Solution Parameters
		params.D =  [1, 1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+ and Cl-.
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
        params.Nseg  =  21; 		% Number of segments between the pipette electrode and tip
		params.V_hold = 100e-3; 	% Initial holding voltage
        params.Rseries = 0e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
        
	% Cell and Membrane parameters
        params.membrane_x = 0.001e-6;   			% Position of membrane patch back from the tip (in m)
		params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_perm = [ 0 0 ]'; 		% Membrane permability to each ion in m/s
        
    % Plot response for 3 pipettes with R = 5 MOhms.
    figure(1); clf;

        t_val = 10.^[-9:0.1:2] ; % Time points
        idx_V = find(t_val < 50e-9);
        idx_rho = find(t_val < 50e-3);
        sf = 1.5; % Scale each pipette diameter by a factor of sf
        r_tip = [ 1/sf, sf, 1] * params.r_tip; % Tip Diameters
        theta = [ sf, 1/sf, 1] * 10 /180 *pi;% Tip Angles
        col = 'gbk'; % Color for each pipette tip size
       
        % Treat each tip
        for j=1:length(r_tip)
            
            % Set pipette tip radius and angle
            paramst = params;
            paramst.r_tip = r_tip(j); % Set Tip Radius
            paramst.theta = theta(j); % Set Tip Angle
            paramst = initialize_system(paramst);
            
            % Calculate scalings
            ltip(j) = r_tip(j) / sin(theta(j)/2);     
            tau_rho(j) = ltip(j)^2/params.D(1);
            R(j) = sum(resistance(paramst));
            
            % Equilibrate initial state
            [paramst, results] = evolve_backward_euler(paramst, inf);

            % Set membrane permeability for IK+ * R = 5mV
             % Calculate membrane permeability for an initial current of 2nA K+
            I_mem = [1e-9, 0 ]'; % Initial current of 10nA by K+
            P_GHK = calculate_PGHK_from_I(I_mem, paramst.rho_pip, paramst.rho_bath, paramst.V_hold, paramst.z, paramst.A(paramst.membrane_N), paramst.T);
            
            paramst.membrane_perm = P_GHK;
            
            % Do time evolution
            [paramst, results] = evolve_backward_euler(paramst, t_val);
            
            % And find steady state
            [params_f, results_f] = evolve_backward_euler(paramst, inf);

            
            % Find time for density change
            rho = mean(results.rho_intra,1);
            rho_min = mean(results_f.rho_intra,1);
            rho_max = mean(results.rho_intra(:,1),1);
            rho_mid = ( rho_max + rho_min * (exp(1)-1) ) / exp(1);
            [rho_asc, idx_asc] = unique(rho);
            t_mid_rho(j) = interp1(rho_asc, t_val(idx_asc), rho_mid);
            
            I = results.I_pip;
            I_min = results_f.I_pip;
            I_max = max(results.I_pip);
            
            % Plot
            subplot(2,2,1); plot(t_val(idx_rho)*1e3, results.I_pip(idx_rho) * 1e9, col(j)); hold on; % I versus t
            subplot(2,2,2); semilogx(t_val/tau_rho(j), (I - I_min) / (I_max - I_min), col(j)); hold on; % I versus tscaled   
            
            subplot(2,2,3); plot(t_val(idx_rho)*1e3, rho(idx_rho), col(j)); hold on;  % rho versus t
            subplot(2,2,4); semilogx(t_val/tau_rho(j), (rho - rho_min)/ (rho_max - rho_min), col(j)); hold on; % rho versus t_scaled
            
            fprintf(1, '\n d_tip = %.2f microns, theta = %.2f degrees, R = %.2f MOhms, P_K = %.2f * 10^{-5} m/s, tau_rho %.2f msec, tmid_rho = %.2f msec', 2*r_tip(j)*1e6, theta(j)*180/pi, R(j)/1e6, paramst.membrane_perm(1)*1e5, tau_rho(j) * 1e3, t_mid_rho(j)  * 1e3);
       
        end

        
        % Label Plots
        subplot(2,2,1);
                xlabel('t (ms)');
                ylabel('I_{pip}(t) (nA)');
                set(gca,'XTick',[0 20 40 ]);
                set(gca,'YTick',[0.8 0.85 0.9 0.95]);
                axis([0 40 0.79 0.97])
                set(gca,'TickDir','out'); 
                box off
                    
       subplot(2,2,2);
                % axis([ 10^(-1.5), 10^(1.5), 0 1]);
                axis([ 1e-2 1e2 0 1])
                xlabel('t / \tau_{diff}');
                ylabel('\DeltaI_{pip}(t) / \DeltaI_{pip}(\infty)) ');
                set(gca,'XTick',[1e-2 1e-0  1e2 ]);
                set(gca,'YTick',[0, 0.5, 1]);
                set(gca,'TickDir','out');
                box off
                
        subplot(2,2,3);
                %axis([ 10^(-5), 10^(1), 125 150]);
                axis([0 40 127 152])
                xlabel('t (ms)');
                ylabel('\rho(t) (mM)');
                set(gca,'XTick',[0 20 40 ]);
                set(gca,'YTick',[130 140 150]);
                set(gca,'TickDir','out');
                box off;
                
        subplot(2,2,4);
                %axis([ 10^(-3), 10^(3), 0 1]);
                axis([1e-2 1e2 0 1])
                xlabel('t / \tau_{dif} ');
                ylabel(' \Delta\rho(t) / \Delta\rho(\infty))');
                set(gca,'XTick',[1e-2 1e-0  1e2 ]);
                set(gca,'YTick',[0 0.5 1]);
                set(gca,'TickDir','out');
                box off;
                
        % Save Figure
        save_figure('Fig_Tau_Pipette_1B',[12 14]);
                
    % Plot Resistance versus 
    % Tip Diameter
        d_min = 2e-7; d_max = 1e-5; Nd = 200;
        d_tip = d_min * 10.^ ( log10(d_max/d_min) * [0:Nd]/Nd   ) ; % 
        
        theta_min = 5; theta_max = 20; Nt = 200;
        theta = theta_min * 10.^( (log10(theta_max/theta_min)) * [0:Nt]'/Nt);
    
        % Square Versions
        d_tip_s = ones(size(theta)) * d_tip;
        theta_s = theta * ones(size(d_tip));
    
        % Resistance
        paramst = initialize_system(params);
        sigma = (paramst.F^2 / (paramst.R * paramst.T) ) * sum( paramst.z .* paramst.z .* paramst.rho_pip .* paramst.D);
        R  = 1./ ( sigma * pi * d_tip_s .* tand(theta_s/4));
        
        % Tau_diffusion
        tau_dif = ( (d_tip_s/2) ./ sind (theta_s/2) ).^2 / params.D(1);
                
        
        figure(2); clf;
        
         subplot(1,2,1); % Resistance 
              imagesc(log10(R));
              axis image;
              set(gca,'YDir','normal');
              xlabel('d_{tip} (\mu m)');
              ylabel('\theta (\circ)');
              xtick = [0.2 0.5 1 2 5 10];
              xtick_pos = interp1(d_tip , 0.5 + [0:Nd],xtick * 1e-6,'linear','extrap');
              set(gca,'XTickLabel',xtick, 'XTick', xtick_pos);
              ytick = [5 7 10 15 20];
              ytick_pos = interp1(theta , 0.5 + [0:Nd],ytick );
              set(gca,'YTickLabel',ytick, 'YTick', ytick_pos);
              set(gca,'TickDir','out');
              colormap(jet(256));
              
              h = colorbar('SouthOutside');
              ctick = [0.2 0.5 1 2 5 10 20];
              set(h,'XTick', 6 + log10(ctick), 'XTickLabel', ctick,'TickDir','out');
              title('R (M\Omega)')          
              
        subplot(1,2,2); % Tau_diff
                imagesc(log10(tau_dif))
        
                axis image;
                set(gca,'YDir','normal');
                xlabel('d_{tip} (\mu m)');
                ylabel('\theta (\circ)');
                xtick = [0.2 0.5 1 2 5 10];
                xtick_pos = interp1(d_tip , 0.5 + [0:Nd],xtick * 1e-6,'linear','extrap');
                set(gca,'XTickLabel',xtick, 'XTick', xtick_pos);
                ytick = [5 7 10 15 20];
                ytick_pos = interp1(theta , 0.5 + [0:Nd],ytick );
                set(gca,'YTickLabel',ytick, 'YTick', ytick_pos);
                set(gca,'TickDir','out');
              
                h = colorbar('SouthOutside');
                ctick = [ 1e-4 1e-3 1e-2 1e-1 1];
                    set(h,'XTick', log10(ctick),'TickDir','out');
                colormap(jet(256));
                title('\tau_{diff} (s)')
                
        save_figure('Fig_Tau_Pipette_2B',[10 8]);     