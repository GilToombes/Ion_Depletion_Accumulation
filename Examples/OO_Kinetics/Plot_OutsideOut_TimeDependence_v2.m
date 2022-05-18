%
% function Example_OutsideOut_TimeDependence
%
%	Example of a calculation for the inside-out configuration.
%       K-channel time dependent pulse.
% 

function Figure_OutsideOut_TimeDependence

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
    
	% Initialize all remaining parameters
    params = initialize_system(params);
    
    % Calculate membrane permeability for an initial current of 2nA K+
    I_mem = [1e-9, 0 ]'; % Initial current of 10nA by K+
    P_GHK = calculate_PGHK_from_I(I_mem, params.rho_pip, params.rho_bath, params.V_hold, params.z, params.A(params.membrane_N), params.T);
    
    % Report Access Resistance and Series Resistance Correction
   	fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', sum(resistance(params)/1e6), params.Rseries/1e6)
    fprintf('\n Membrane Permeability = [%e, %e] m/s', P_GHK);
    
    
	% Equilibrate initial state
	[params, results] = evolve_backward_euler(params, inf);
        
    % Need to catch each time step to plot relevant quantities
    params.membrane_perm = P_GHK; %[ 9.06e-5, 0]' / 1.645 ; % Make the membrane K-selective
    %t_val = 1e-9 * 10.^[-1:0.2:11];       % Time Steps
    t_range1 = [1:10,11:10:501] * 1e-9; % 
    idx1 = [1:length(t_range1)];
    t_scale1 = [0, 500];
    t_range2 = [0.6e-3, 1:10:99, 100:25:1000, 1000:500:5000] * 1e-6;
    idx2 = length(t_range1) + [1:length(t_range2)] ;
    t_scale2 = [0.6/1000, 5];
    t_range3 = [6:10, 20:10:95, 100:50:1000] * 1e-3;
    idx3 = length(t_range1) + length(t_range2) + [1:length(t_range3)];
    t_scale3 = [5, 200];
    
    t_val = [t_range1, t_range2, t_range3];    
    
    x_val = params.x([1:(params.Nseg-1)]);
    x_tip = params.r_tip / tan( params.theta/2) ;
    x_val_J = 1./ (0.5 ./ (x_val-x_tip) + 0.5 ./(params.x([2:params.Nseg]) - x_tip)) + x_tip;
    Ne = length(t_val);
    Nl = size(params.J,2); % Number of points inside pipette
    idxl = [1:Nl];
    JK = zeros(Ne, Nl);
    JCl =zeros(Ne, Nl);
    rho_K = zeros(Ne, Nl);
    rho_Cl = zeros(Ne, Nl);
    V = zeros(Ne, Nl);
    
    for j = 1:length(t_val)
        [params, results] = evolve_backward_euler(params, t_val(j));
        JK(j,:) = params.J(1,:);
        JCl(j,:) = params.J(2,:);
        rho_K(j,:) = params.rho(1,idxl);
        rho_Cl(j,:)= params.rho(2,idxl);
        V(j,:) = params.V(idxl);
    end
    
    Nmem = params.Nseg-1;
    Imem = JK(:,Nmem) *params.F*1e9;
    Ielec = (JK(:,1)-JCl(:,1)) *params.F*1e9;
    
    
    
    figure(1); clf
    % Plot I, V, rho, delta_mu versus time.
    % Split timebase into 3 parts
    
    fig_dim = [10 8.3];
    save_figure('Fig_Example_OutsideOut_TimeDependence_A',fig_dim)
    
    % Plot Positions
        xp = [ 0.1 0.40 0.70 ]; % Position of 3 columns
        yp = [0.08, 0.32, 0.57, 0.75]; % Position of rows
        yh = [ 0.2, 0.2,  0.15, 0.2]; % Height of rows
        xw = [ 0.25, 0.25, 0.25]; % Width of columns
        
    %       subplot(1,3,3); hold on; plot(x_val*1e6,  log(rho_K(idx(j),:)/150) + (V(idx(j),:)-1e-1) / (params.R * params.T / params.F ));            axis([-20 0 -0.4 0.4])
    %       plot(x_val*1e6,  log(rho_Cl(idx(j),:)/150) - (V(idx(j),:)-1e-1) / (params.R * params.T / params.F ));            
       
     
        % Currents on Top
        subplot('position', [xp(1) yp(4) xw(1), yh(4)]);
        y_range = [0 1];
        plot(t_range1*1e9, Ielec(idx1) ,'k--'); hold on;      plot(t_range1*1e9,Imem(idx1),'k'); 
        axis([t_scale1, y_range]);
        box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out';
        ylabel('I (nA)'); 
        
        subplot('position', [xp(2) yp(4) xw(2), yh(4)]);
        plot(t_range2*1e3, Ielec(idx2) ,'k--'); hold on;      plot(t_range2*1e3,Imem(idx2),'k'); 
        axis([t_scale2, y_range]);
        box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
        
        subplot('position', [xp(3) yp(4) xw(2), yh(4)]);
        plot(t_range3*1e3, Ielec(idx3) ,'k--'); hold on;      plot(t_range3*1e3,Imem(idx3),'k'); 
        box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
        axis([t_scale3, y_range]);
        
        % Voltages next
        subplot('position', [xp(1) yp(3) xw(1), yh(3)]);
        y_range = [94 100];
        plot(t_range1*1e9, V(idx1,Nmem)*1e3 ,'k'); 
        axis([t_scale1, y_range]);
        box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out'; 
        ylabel('V (mV)'); 
        
        subplot('position', [xp(2) yp(3) xw(2), yh(3)]);
        plot(t_range2*1e3, V(idx2, Nmem)*1e3 ,'k');
        axis([t_scale2, y_range]);
        box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
        
        subplot('position', [xp(3) yp(3) xw(2), yh(3)]);
        plot(t_range3*1e3, V(idx3,Nmem)*1e3 ,'k');
        box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
        axis([t_scale3, y_range]);
        
        % Concentrations third
        subplot('position', [xp(1) yp(2) xw(1), yh(2)]);
        y_range = [125 151];
        plot(t_range1*1e9, rho_K(idx1,Nmem) ,'b'); hold on;      plot(t_range1*1e9,rho_Cl(idx1,Nmem),'g--'); 
        axis([t_scale1, y_range]);
        box off; h = gca;  h.XAxis.Visible = 'off';     h.TickDir='out';
        ylabel('\rho (mM)'); 
        
        subplot('position', [xp(2) yp(2) xw(2), yh(2)]);
        plot(t_range2*1e3, rho_K(idx2,Nmem) ,'b');hold on;      plot(t_range2*1e3,rho_Cl(idx2,Nmem),'g--'); 
        axis([t_scale2, y_range]);
        box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
        
        subplot('position', [xp(3) yp(2) xw(2), yh(2)]);
        plot(t_range3*1e3, rho_K(idx3,Nmem) ,'b'); hold on;      plot(t_range3*1e3,rho_Cl(idx3,Nmem),'g--'); 
        box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
        axis([t_scale3, y_range]);
        
        % Chemical Potentials fourth
        subplot('position', [xp(1) yp(1) xw(1), yh(2)]);
        mu_K = log(rho_K(:,Nmem)/params.rho_pip(1)) + (V(:,Nmem)-params.V_hold) / (params.R * params.T / params.F );  
        mu_Cl = log(rho_Cl(:,Nmem)/params.rho_pip(2)) -(V(:,Nmem)-params.V_hold) / (params.R * params.T / params.F );  
        y_range = [-0.4 0.4];

        plot(t_range1*1e9, mu_K(idx1) ,'b'); hold on;      plot(t_range1*1e9,mu_Cl(idx1),'g--'); 
        axis([t_scale1, y_range]);
        box off; h = gca;  h.TickDir = 'out'; %h.XAxis.Visible = 'off';    
        ylabel('\Delta \mu (k_B T)'); 
        set(gca,'XTick',[250 500]);
        
        subplot('position', [xp(2) yp(1) xw(2), yh(2)]);
        plot(t_range2*1e3, mu_K(idx2) ,'b');hold on;      plot(t_range2*1e3,mu_Cl(idx2),'g--'); 
        axis([t_scale2, y_range]);
        box off ; h = gca;  h.YAxis.Visible = 'off'; h.TickDir = 'out';  % h.XAxis.Visible = 'off';    
        set(gca,'XTick',[2.5 5]);
        
        subplot('position', [xp(3) yp(1) xw(2), yh(2)]);
        plot(t_range3*1e3, mu_K(idx3) ,'b'); hold on;      plot(t_range3*1e3,mu_Cl(idx3),'g--'); 
        box off ; h = gca; h.YAxis.Visible = 'off';  % h.XAxis.Visible = 'off';    
        axis([t_scale3, y_range]);
        set(gca,'TickDir','out');
        
        
    save_figure('Fig_Example_OutsideOut_TimeDependence_A',fig_dim)
        
        
    figure(2);  clf;
    
    fig_dim = [10 8.5];
    save_figure('Fig_Example_OutsideOut_TimeDependence_B',fig_dim)
    
    % Plot pipette state at 3 timepoints.
    idx_time = [length(t_range1) , length(t_range1) + length(t_range2), length(t_val) ] ; % Timepoints to plot

   % Define plots 
        xp = [ 0.1 0.42 0.74 ]; % Position of 3 columns
        yp = [0.08, 0.28, 0.48, 0.68, 0.9 ]; % Position of rows
        yh = [ 0.18, 0.18,  0.14, 0.18, 0.09]; % Height of rows
        xw = [ 0.25, 0.25, 0.25]; % Width of columns
    
    for j = 1:length(idx_time);

    x_range = [-20 0]; % Range from -20 microns to 0

    % Plot Pipette
    subplot('position',[ xp(j), yp(5), xw(j), yh(5)]);
        ytip = params.r_tip*1e6;
        yback = ytip - x_range(1) * tan(params.theta/2);
        plot(x_range, [yback, ytip],'k','Linewidth',1); hold on
        plot(x_range, -1 * [yback, ytip],'k','Linewidth',1); hold on
        
%        for j=1:length(x_val)
%            plot( x_val(j) * 1e6 + [0 0], [-1 1] * (ytip - x_val(j)*1e6 * tan(params.theta/2) ),'k:'); hold on
%        end
        axis([x_range, -2.1 2.1]) ;axis image;  axis off;
        
    % Plot Current
    subplot('position',[ xp(j), yp(4), xw(j), yh(4)]);
    plot(x_val*1e6,  (JK(idx_time(j),:) - JCl(idx_time(j),:)) * params.F * 1e9,'k'); hold on;      
    plot(x_val*1e6,  JK(idx_time(j),:) * params.F * 1e9,'b'); hold on;  
        plot(x_val*1e6,  -1 * JCl(idx_time(j),:) * params.F * 1e9,'g--'); 
        axis([x_range 0 1]);
        box off; h = gca;  h.XAxis.Visible = 'off';     h.TickDir='out';
        if (j==1) ylabel('I (nA)'); else ; h.YAxis.Visible = 'off'; end
        
    % Plot V
    subplot('position',[ xp(j), yp(3), xw(j), yh(3)]);
        plot(x_val*1e6,  V(idx_time(j),:)*1e3,'k'); 
        axis([x_range 95 100]);
        box off; h = gca;  h.XAxis.Visible = 'off';     h.TickDir='out';
        if (j==1) ylabel('V (mV)'); else ; h.YAxis.Visible = 'off'; end
        
    % Plot Rho
    subplot('position',[ xp(j), yp(2), xw(j), yh(2)]);
    plot(x_val*1e6, rho_K(idx_time(j),:),'b'); hold on;
    plot(x_val*1e6, rho_Cl(idx_time(j),:),'g--'); hold on;
    axis([x_range 125 150]);
    box off; h = gca;  h.XAxis.Visible = 'off';     h.TickDir='out';
    if (j==1) ylabel('\rho (mM)'); else ; h.YAxis.Visible = 'off'; end
    
    % Plot Mu
    subplot('position',[ xp(j), yp(1), xw(j), yh(1)]);
    mu_K =  log(rho_K(idx_time(j),:)/ params.rho_bath(1)) +  (V(idx_time(j),:)-params.V_hold) / (params.R * params.T / params.F );  
    mu_Cl = log(rho_Cl(idx_time(j),:)/ params.rho_bath(2)) - (V(idx_time(j),:)-params.V_hold) / (params.R * params.T / params.F );  
    
    plot(x_val*1e6, mu_K,'b'); hold on;
    plot(x_val*1e6, mu_Cl,'g--'); hold on;
    axis([x_range -0.4 0.24]);
    box off; h = gca;  h.TickDir='out'; % h.XAxis.Visible = 'off';    
    if (j==1) ylabel('\Delta \mu (k_BT)'); else ; h.YAxis.Visible = 'off'; end
    xlabel('x (\mu m)')  ;  
    
    end
    
    save_figure('Fig_Example_OutsideOut_TimeDependence_B',fig_dim)
    
    
    
   