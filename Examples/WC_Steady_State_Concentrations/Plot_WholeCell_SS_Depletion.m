% function Figure6_WholeCell_SS_Depletion
%
%       Illustrate that ion accumulation/depletion depends on individual
%       ionic currents rather than net current.
%
%       Whole-Cell
%           Pipette - Raccess = 5 MOhms divided into 60 segments; 
%                       d_tip =  0.64 microns ; 	theta_tip = 10 degrees
%           Pipette Solution  = Intra-cellular = 150mM K-gluconate
%           Bath Solution     = Extra-Cellular = 150mM NaCl
%
%       
%       Plot rho,V versus position 'x' from pipette tip
%       for a pipette containing 150mM KCl, 5 MOhms
%       Show   ICl- = 1nA;
%	       IK+  = 1nA;
%              I_Na+ = -1nA;
%	       Ik+ = 1nA; I_Na+ = -1nA
%

function Figure_6_WholeCell_SS_Depletion

	% Configuration parameters
		params.config = 'whole-cell';      % Whole-cell 
        params.Constraints.type = 'basic';  % Option specifying which parameters are variable
                                            %       basic   - voltage, molecular concentrations and currents
	% Solution Parameters
		params.D =  [1, 0.33, 0.682, 1.0382]' * 1.957e-9; 	% Diffusion coefficients of K+, gluc-, Na+ and Cl-.
		params.z =  [1, -1, 1, -1]';	  	% Relative charge of K+, Na+ and Cl-.
		params.rho_bath = [0, 0, 150, 150  ]';  	% 150mM NaCl in bath
		params.rho_pip =  [150, 150, 0, 0   ]';  	% 150mM K-gluc in pipette
        
	% Pipette parameters
        R_pip = 5e6 ; % Pipette resistance of 5 MOhms
		params.r_pip = 0.5e-3; % Pipette inner radius at electrode (in m)
        params.theta = 10 * pi/180; % Tip Angle in radians
        [tip_diameter, sigma] = calculate_tip_diameter(params, R_pip); % Calculate tip radius
        tip_diameter = tip_diameter * 5.23/5 ; % Account for small extra resistance between pipette tip and cell membrane
        params.r_tip = tip_diameter/2;% Tip radius of pipette (in m)
        params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
        params.Nseg  =  61; 		% Number of segments between the pipette electrode and tip
		params.V_hold = 0e-3;       % Initial holding voltage
        params.Rseries = 0e6;       % Series resistance correction in Ohms
        params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
                                    
	% Cell and Membrane parameters
        params.membrane_x = 1e-10;   			% Position of membrane patch back from the tip (in m)
		params.membrane_func = @membrane_current_constant; 	% Function describing membrane permeability to ions
		params.membrane_J = [ 0 0 0]'; 		% Membrane molecular current for each ion in mole/s
        
       
   % Cell and Membrane parameters
     	
        params.A_cell = 1e-9; % Cell Area
        params.A_cell = 0.767e-9; % Cell Area
        params.Lambda_Cell = 2e-15 ; % Cell volume
            
        % Calculate positions and segments
            Npip  = 60 ; % Number of pipette segments
            Ncell = 5  ; % Number of segments from pipette to cell surface
            Nouter = 5;  % Number of segments from outer membrane surface to external pipette.
            
            % Pipette segments
            omega_pip = 4 * pi * sin(params.theta/4)^2; %   Solid angle
            dtip =  params.r_tip / tan(params.theta/2); % Distance from source to tip.
            delec = params.r_pip / tan(params.theta/2); % Distance from source to electrode
            xpip = dtip - 1./(1/delec + (1/dtip - 1/delec) * ([1:Npip]-1)/(Npip-1)); %1/r spacing for pipette
            Apip = omega_pip * (dtip-xpip).^2;
            
            % Cell segments
            Rcell = 3 * params.Lambda_Cell / params.A_cell; % Cell radius
            omega_cell = params.A_cell / Rcell^2;
            rinner = (omega_pip / omega_cell)^0.5 * dtip;
            router = Rcell;
            xcell = [ 1./ (1/rinner + (1/router-1/rinner)*[1:Ncell]/Ncell) ]-rinner; % 1/r spacing for cell
            Acell = omega_cell * (xcell+rinner).^2;
            
            % Bath segments
            rmem = Rcell ;
            rbath = params.r_pip;
            xbath = [ 1./ (1/rmem + (1/rbath - 1/rmem)* [1:Ncell]/Ncell) ] - rinner; % 1/r spacing out to bath pipette
            Abath = omega_cell * (xbath+rinner).^2;
            
            % Combine segments.
            params.x = [xpip, xcell, xbath];
            params.A = [Apip, Acell, Abath];
            params.membrane_N = Npip + Ncell ; % Position of membrane
            params.Nseg  =  length(params.x); 		% Number of segments between the pipette electrode and tip

        
        
	% Initialize all remaining parameters
    params = initialize_system(params);
    
    % Report Access Resistance and Series Resistance Correction
    fprintf('\n Pipette Solution Conductance = %.2f S/m. Tip Diameter = %.2f microns.', sigma, params.r_tip * 2e6);
    fprintf('\n Access Resistance = %.2f MOhms. Series Resistance Correction = %.2f MOhms.', resistance(params)/1e6, params.Rseries/1e6);

    figure(1); clf;
    figure('DefaultAxesFontSize',6);
    pw = 17.5/2;    % Width of Example in cm
    ph = 7.5;         % Hight of Example in cm; 
    
        % IK = 1nA
        Imem = 1e-9;
        params.membrane_J = Imem/params.F * [ 1, 0 , 0, 0]' ; 		% Membrane permability to each ion in m/s
        plot_it2(params,2); 
        save_figure('Figure6_WholeCell_K',[pw, ph]);  
        
        %ICl- = 1nA
        figure(1); clf;
        params.membrane_J = Imem/params.F * [ 0, 0, 0, -1]' ; 		% Membrane permability to each ion in m/s
        plot_it2(params,2); 
        save_figure('Figure6_WholeCell_Cl',[pw, ph]);  
        
        
        % INa+ = -1nA
        figure(1); clf;
        params.membrane_J = Imem/params.F * [ 0, 0, -1, 0]' ; % Membrane permability to each ion in m/s
        plot_it2(params,3);
        save_figure('Figure6_WholeCell_Na',[pw, ph]);  
        
        % INa+ = -1nA, IK+ = 1nA
        figure(1); clf;
        params.membrane_J = Imem/params.F * [ 1, 0, -1, 0]' ; % Membrane permability to each ion in m/s
        plot_it2(params,4);  
        save_figure('Figure6_WholeCell_K_and_Na',[pw, ph]);            
            

    
    function plot_it2(params,pos)
        
       % Inputs :
       %            params -> state of pipette
       %            pos = 1, 2, 3, or 4 -> column to plot in
       
       
       % Positions/Scale of panels for left side
       left_rhoV =  0.13 ; % Left horizontal corner of rho/V cell 
       width_rhoV = 0.25; % Width of 
       xrange_rhoV = [-18 20]; % X-axis range of rho/V cell plots
       
       bottom_rho = 0.03; % Bottom vertical of rho cell plot
       height_rho = 0.35; % Height of rho cell plot
       yrange_rho = [0 180]; % Y-range of rho plot
       
       bottom_V = 0.45; % Bottom of vertical of V cell plot
       height_V = 0.10; % Height of V cell plot
       yrange_V = [-3.5 5.1]; % Y-range of V plot
       
       %bottom_mu = 0.75; % 
       %height_mu = 0.15;
       %yrange_mu = [-0.3 0]; % Y-range of mu plot
       
       % Right side Plots of Concentration and Voltage Change at membrane verus I
       left_Ohm =  0.6 ; % Left horizontal corner of rho/V cell 
       width_Ohm = 0.23; % Width of 
       xrange_I = [0 1.1]; % X-axis range in 1 nA
       
       bottom_rhoOhm = 0.45; % Bottom vertical of rho cell plot
       height_rhoOhm = 0.48; % Height of rho cell plot
       yrange_rhoOhm = [0 180]; % Y-range of rho plot
       
       bottom_VOhm = 0.1; % Bottom of vertical of V cell plot
       height_VOhm = 0.23; % Height of V cell plot
       yrange_VOhm = [-3.5 5.1]; % Y-range of V plot
       
       
       if (nargin > 1)
           switch pos
               case 1
                    yrange_rho = [0 185]; % Y-range of rho plot
                    yrange_V = [-4 0]; % Y-range of V plot
                    yrange_rhoOhm = [0 190]; % Y-range of rho plot
                    yrange_VOhm = [-0.1 5.1]; % Y-range of V plot
            
               case 2
                    yrange_rho = [0 185]; % Y-range of rho plot
                    yrange_V = [-4 0.2]; % Y-range of V plot
                    yrange_rhoOhm = [0 190]; % Y-range of rho plot
                    yrange_VOhm = [-0.1 5.1]; % Y-range of V plot
                    
               case 3
                    yrange_rho = [0 185]; % Y-range of rho plot
                    yrange_V = [-0.2 5.1]; % Y-range of V plot
                    yrange_rhoOhm = [0 190]; % Y-range of rho plot
                    yrange_VOhm = [-5.1 0.1]; % Y-range of V plot
    
              case 4
                    yrange_rho = [0 185]; % Y-range of rho plot
                    yrange_V = [-0.2 5.1]; % Y-range of V plot
                    yrange_rhoOhm = [0 190]; % Y-range of rho plot
                    yrange_VOhm = [-5.1 0.1]; % Y-range of V plot
    
                    
                    
           end
       end
       
       
       Raccess = max(resistance(params)); % Calculate Access Resistance
       [params, results] = evolve_backward_euler(params, inf);        % Solve for steady-state concentraitons
       I = params.J(:,params.membrane_N) * params.F * 1e9 .* params.z; % Check Membrane Current
       fprintf(1,'\n IK+ = %.3f nA, Igluc- = %.2f nA, INa+ = %.2f nA, ICl- = %.2f nA', I)

      
     
       
       % Rescale distance from pipette tip
       xval = 1e6 * params.x; % Xposition in microns

       % Rescale pipette distances
       idx_pip = [1:60]; % Pipette indices
       xpip = xval(idx_pip); % Pipette x values
       lpip = 18; % Pick scaling length in microns
       xpiprs = lpip * ( lpip./(lpip - xpip) - 1);
       xval(idx_pip) = xpiprs;

       % Rescale cell distances
       idx_cell = [61:65]; % Cell Indices
       idx_pc = [idx_pip, idx_cell];
       lcell = 14; % Pick cell length scale 
       xcell = xval(idx_cell);
       xval(idx_cell) = xcell/max(xcell) * lcell;

       % Rescale bath distances
       idx_bath = [66:70];
       lbath = 6; % Bath length scale
       xval(idx_bath) = lcell + (idx_bath- min(idx_bath))/(max(idx_bath) - min(idx_bath)) * lbath;
       
       
       %Concentrations versus position     
       subplot('Position',[ left_rhoV, bottom_rho, width_rhoV, height_rho]); 
            plot(xval(idx_pc), params.rho(1,idx_pc), 'b'); hold on
            plot(xval(idx_bath), params.rho(1,idx_bath), 'b'); hold on
            plot(xval(idx_pc), params.rho(2,idx_pc), 'Color',[0.8, 0, 0]); 
            plot(xval(idx_bath), params.rho(2,idx_bath), 'Color',[0.8, 0, 0]); 
            plot(xval(idx_pc), params.rho(3,idx_pc), 'Color', [0.8, 0.8, 0], 'Linewidth',1); 
            plot(xval(idx_bath), params.rho(3,idx_bath), 'Color', [0.8, 0.8, 0], 'Linewidth',1); 
            plot(xval(idx_pc), params.rho(4,idx_pc), 'Color', [0.6, 0.6, 0.6], 'Linewidth',1,'LineStyle','--'); 
            plot(xval(idx_bath), params.rho(4,idx_bath), 'Color', [0.6, 0.6, 0.6], 'Linewidth',1,'LineStyle','--'); 
            axis([xrange_rhoV, yrange_rho])
            set(gca,'XTick',[]);
            set(gca,'YTick',[0 30 60 90 120 150 180]);
            ylabel('\rho (mM)');
            set(gca,'TickDir','out'); box off

       % Voltage versus position
       subplot('Position',[ left_rhoV, bottom_V, width_rhoV, height_V]);
            V = params.V - params.V_hold  ;
            plot(xval(idx_pc), V(idx_pc)*1e3, 'k' ); hold on
            plot(xval(idx_bath), V(idx_bath)*1e3, 'k' ); hold on
            axis([xrange_rhoV, yrange_V])
            set(gca,'XTick',[]);
            ylabel('V (mV)');
            set(gca,'TickDir','out'); box off;
    
       %  Mu versus position
%       subplot('Position',[ left_rhoV, bottom_mu, width_rhoV, height_mu]);
%            mu_K =    log(params.rho(1,:)/params.rho(1,1)) + (params.V - params.V_hold) * params.F / (params.R * params.T);
%            mu_gluc = log(params.rho(2,:)/params.rho(2,1)) - (params.V - params.V_hold) * params.F / (params.R * params.T);
%            plot(xval(idx_pc),   mu_K(idx_pc), 'b' ); hold on
%            plot(xval(idx_pc),   mu_gluc(idx_pc), 'Color', [0.8 0 0] ); hold on
%            axis([xrange_rhoV, yrange_mu]);
%            set(gca,'XTick',[]);
%            ylabel('\Delta\mu (kT)');
%            set(gca,'TickDir','out'); box off;
    
      
            mem_J = params.membrane_J; % Maximum membrane current
            sf = [0.25, 0.5, 0.75, 1];   % Fractional current levels to calculate
            Np = length(sf)+1;
            rho_numeric = params.rho_pip * ones(1,Np);
            rho_analytic = rho_numeric;
            V_numeric = zeros(Np,1);
            V_analytic = zeros(Np,1);
            I = zeros(Np,1);
            
            % Calculate delta_rho for each current level
            for j=1:length(sf)
        
                params.membrane_J = mem_J * sf(j) ; 		% Membrane current 
                [params, results] = evolve_backward_euler(params, inf); % SS concentrations
                Ia = params.J(:,params.membrane_N) * params.F .* params.z; % Current
                sigma = sum( params.rho_pip .* params.D .* params.z.^2  .* params.F^2 / (params.R * params.T) );
                rhoz = params.rho_pip .* params.z .* params.z;
                rho_analytic(:,j+1) = params.rho_pip - 1 * Raccess * sigma ./ (params.F * params.z .* params.D) .*  (Ia - rhoz .* params.D * sum( Ia ./ params.D) / sum(rhoz)   );
                rho_numeric(:,j+1) = results.rho_intra  ;
                V_analytic(j+1) = Raccess * sum( Ia ./ params.D) * sum(rhoz .* params.D) / sum(rhoz)  *1e3 ;
                V_numeric(j+1) = results.V_mem * -1e3; 
                V_ohm(j+1) = sum(Ia) * Raccess *1e3;
                I(j+1) = max(abs(Ia)) * 1e9;
                
         %       D_eff = sum(params.D .* rhoz) / sum(rhoz);
         %       rho_eff = sum(rhoz)/2;
         %       VT = params.R * params.T / params.F;
         %      Va = Ia * Raccess ./ (params.D / D_eff);
         %       V_analytic(j+1) = sum(Va);
         %       rho_analytic(:,j+1) = params.rho_pip -  (2 * rho_eff * Va./ (params.z * VT) - (params.z .* params.rho_pip/VT) * V_analytic(j+1) ); 
         %       V_analytic(j+1) = V_analytic(j+1) * 1e3;
            end

            subplot('Position',[ left_Ohm, bottom_rhoOhm, width_Ohm, height_rhoOhm]); 
            plot( I, rho_analytic(1,:),'b','Marker','*','LineStyle','--'); hold on
            plot( I, rho_analytic(2,:),'Color',[0.8, 0, 0],'Marker','*','LineStyle','--'); hold on
            plot( I, rho_analytic(3,:),'Color',[0.8, 0.8, 0],'Marker','*','LineStyle','--'); hold on
            plot( I, rho_analytic(4,:),'Color',[0.6, 0.6, 0.6],'Marker','*','LineStyle','--'); hold on
            plot( I, rho_numeric(1,:),'b','Marker','s'); hold on
            plot( I, rho_numeric(2,:),'Color',[0.8, 0, 0],'Marker','s'); hold on
            plot( I, rho_numeric(3,:),'Color',[0.8, 0.8, 0],'Marker','s'); hold on
            plot( I, rho_numeric(4,:),'Color',[0.6, 0.6, 0.6],'Marker','s'); hold on

            axis([xrange_I, yrange_rhoOhm]);
            set(gca,'XTick',[0 0.25 0.5 0.75 1]);
            set(gca,'TickDir','out'); box off;
            ylabel('\rho_{mem} (mM)');
            set(gca,'YTick',[0 30 60 90 120 150 180]);
            
            subplot('Position',[ left_Ohm, bottom_VOhm, width_Ohm, height_VOhm]); 
            plot( I, V_analytic,'k','Marker','*','LineStyle','--'); hold on
            plot( I, V_numeric,'k','Marker','s'); hold on
            plot( I, V_ohm,'Color',[0.5 0.5 0.5]);
            set(gca,'XTick',[0 0.25 0.5 0.75 1]);
            set(gca,'TickDir','out'); box off;
            xlabel('|I| (nA)'); 
            ylabel('\DeltaV (mV)');
            set(gca,'TickDir','out'); box off;
            axis([xrange_I, yrange_VOhm]);
            ax = gca; ax.XAxisLocation = 'origin';
            fprintf(1,'\n Analytic rhoK = %.2f mM, rho_gluc = %.2f mM, rhoNa = %.2f mM, rhoCl = %.2f mM',  rho_analytic(:,Np));
            fprintf(1,'\n Numeric  rhoK = %.2f mM, rho_gluc = %.2f mM, rhoNa = %.2f mM, rhoCl = %.2f mM',rho_numeric(:,Np));
            fprintf(1,'\n Numeric  DeltaV = %.2f mV, Analytic = %.2f mV,\n', V_numeric(Np), V_analytic(Np));
        
       
        