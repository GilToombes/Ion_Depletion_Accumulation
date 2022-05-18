% function Figure_Buffering_Calcium_V2
%
%   Illustrate effects of calcium buffering
%           
%           Panel A-C
%               Outside-out patch 
%                       R = 5 MOhms containing 150mM K+
%                       A = 10mM BAPTA + 110mM Cl-, 1pA inward Ca++ current
%                       B = 10mM EGTA + 130mM Cl-,  1pA inward Ca++ current
%                       C = 10mM BAPTA + 110mM Cl-, 110 pA inward Ca++ current
%                       bath - 150mM KCl + 2mM Ca++
%
%               Plot
%                       I versus x (position from membrane)
%                           free calcium; Ca++:buffer 
%                       Concentrations versus x
%                           Free Buffer - blue
%                           Calcium-bound buffer - green
%                           Free calcium -  black
%                       Binding rate versus x
%           
%           Panel D
%                   Free calcium versus |ICa++| for inward current
%                       Black : R = 5 MOhms oo patch no buffer
%                       Blue :  R = 5 MOhms oo patch 10mM EGTA
%                       Green : R = 5 MOhms oo patch 10mM BAPTA
%                       Cyan :  R = 5 MOhms, whole-cell, 10mM EGTA
%                       Red :   R = 5 MOhms, whole-cell, 10mM BAPTA
%                       
%
%       Access Resistance = 5.00 MOhms.
%       Unbuffered solution gives 74.79 uM per pA of Ca++ current.
%
%       10mM BAPTA tau_free = 0.22 microseconds, lambda_free = 13.27 nm
%       10mM EGTA  tau_free = 37.04 microseconds, lambda_free = 171.29 nm
%
%       1.12 microM free Ca++ per picoamp per micron^2 with 10mM EGTA
%       86.79 nM free Ca++ per picoamp per micron^2 with 10mM BAPTA
%
%       Maximum current from free buffer diffusion 82.58 pA.


function Figure_Buffering_Calcium_V2

        params_oo_EGTA  = initialize_state('outside-out', 'EGTA'); % Outside-Out patch 10mM EGTA 
        params_oo_BAPTA = initialize_state('outside-out', 'BAPTA');% Outside-Out patch 10mM BAPTA 
        params_wc_BAPTA = initialize_state('whole-cell',  'BAPTA');% Whole-cell 10mM BAPTA 
        params_wc_EGTA = initialize_state('whole-cell',  'EGTA');  % Whole-cell 10mM EGTA 
        %params_io_BAPTA = initialize_state('inside-out', 'BAPTA');

        % Create unbuffered oo config;
        params_oo_NOBUFFER = params_oo_EGTA;
        for j=1:(params_oo_NOBUFFER.Nseg-1)
            rho_nobuf = [0 0 0 150 150]';
            params_oo_NOBUFFER.rho_pip = rho_nobuf;
            params_oo_NOBUFFER.rho(:,j) = rho_nobuf ; % 150mM KCl in pipette
        end
               
        % Report Access Resistance and Series Resistance Correction
        fprintf('\n Access Resistance = %.2f MOhms.', sum(resistance(params_oo_NOBUFFER)/1e6))
        
        % Calculate Estimates from paper
        % 1. Unbuffered concentration change per picoamp
        params = params_oo_NOBUFFER;
        R_pip = 5e6;
        sigma = sum(params.rho_pip .* params.D .* params.z .* params.z) * params.F^2 / (params.R * params.T); 
        
        fprintf(1, '\n Unbuffered solution gives %.2f uM per pA of Ca++ current.\n',  R_pip * sigma / (params.D(1) * params.z(1) * params.F) * 1e3/1e12);
        
        % 2. Free Ca++ lifetime and distance
        rho_buf = params_oo_BAPTA.rho_pip(2); % Concentration of free buffer
        Kon = [params_oo_BAPTA.Kon, params_oo_EGTA.Kon];
        D = params_oo_BAPTA.D(1); % Diffusion coefficient of Ca++
        tau_free = 1./ (Kon * rho_buf);
        lambda_free = (D ./ (Kon * rho_buf)).^0.5;
        fprintf(1,'\n 10mM BAPTA tau_free = %.2f microseconds, lambda_free = %.2f nm', tau_free(1)*1e6, lambda_free(1) * 1e9); 
        fprintf(1,'\n 10mM EGTA  tau_free = %.2f microseconds, lambda_free = %.2f nm', tau_free(2)*1e6, lambda_free(2) * 1e9); 
        
        % 3. Estimate of delta_rho_Ca++ due to buffer kinetics
        IpA = 1 ; % 1 pA/square micron current density
        drho = IpA ./ ( params.z(1) * params.F * (params.D(1) * Kon * rho_buf).^0.5 ) ;
        fprintf('\n\n %.2f microM free Ca++ per picoamp per micron^2 with 10mM EGTA', drho(2) * 1e3);
        fprintf('\n %.2f nM free Ca++ per picoamp per micron^2 with 10mM BAPTA', drho(1) * 1e6);
        
        %4. Estimate of maximum current 
        Imax = abs(params.z(1)) * params.D(2) * params.F * rho_buf / (R_pip * sigma);
        fprintf('\n\n Maximum current from free buffer diffusion %.2f pA.', Imax * 1e12);
        
        % Plot panel position parameters
        l1 = 0.1; w = 0.35; % Left Column x and width
        l2 = 0.6; % Right Column x
        x_range = [1e-3, 20]; % Range of plot
        y1 = 0.375; y2 = 0.15; y3 = 0.05;
        h1 = 0.09;  h2 = 0.22;  h3 = 0.09;
        yoff = 0.5; % height offset
        
        % Panel 1 = BAPTA-OO-1pA SS
                                    %   x, y, width, height
        plot_options(1).panels = [  l1, yoff+y1, w,  h1; % I versus x
                                    l1, yoff+y2,  w, h2; % C versus x
                                    l1, yoff+y3,  w, h3;]; % d rho_Ca-buff dt
        
        plot_options(1).axes  = [   x_range, 0,      1.05;         % I versus x
                                    x_range, 1e-9,   2e-2 ; % C versus x
                                    x_range, 1e-3,     10;];    % d rho_Ca-buff dt
        
        % Panel 2 = EGTA-OO-1pA SS
                                    %   x, y, width, height
        plot_options(2).panels = [ l2,  yoff+y1, w,  h1; % I versus x
                                   l2,  yoff+y2,  w, h2; % C versus x
                                    l2, yoff+y3,  w, h3;]; % d rho_Ca-buff dt
        
        plot_options(2).axes  = [ x_range, 0,      1.05;         % I versus x
                                  x_range, 1e-9,   2e-2 ; % C versus x
                                  x_range, 1e-3,      10;];    % d rho_Ca-buff dt

        % Panel 3 = BAPTA-OO-110pA SS
                                    %   x, y, width, height
        plot_options(3).panels = [ l1, y1,  w, h1; % I versus x
                                   l1, y2,  w, h2; % C versus x
                                   l1, y3,  w, h3;];% d rho_Ca-buff dt
        
        plot_options(3).axes  = [ x_range, 0,      120;         % I versus x
                                  x_range, 1e-9,   2e-2 ; % C versus x
                                  x_range, 1e-3,      10;];    % d rho_Ca-buff dt

        % Panel 4 - SS rho/I plots
        plot_options(4).panels = [l2, 0.07, w, 0.4]; % Plot position
        
        % Initialiaze figure
        figure(1); clf;

        % Panel 1 : 10mM BAPTA; -1pA Ca++
        plot_rho_vs_x_ss(params_oo_BAPTA,  plot_options(1));  
        
        % Panel 2 : 10mM EGTA ; -1pA Ca++
        plot_rho_vs_x_ss(params_oo_EGTA, plot_options(2));   
        
        % Panel 3 : 10mM BAPTA ; -110pA Ca++
        paramst = params_oo_BAPTA; paramst.membrane_perm = params_oo_BAPTA.membrane_perm * 103 * 1.107 ;
        plot_rho_vs_x_ss(paramst, plot_options(3));  
        
        % Panel 4 : Plot of rho_Ca++ versus current for unbuffered, BAPTA/EGTA oo, BAPTA/EGTA whole-cell
        subplot('Position',plot_options(4).panels(1,:));
                scale_factors = [0.1,0.3, 1,3,10,30,45,65,80,87.5,92, 93.5, 94.25, 95, 95.625, 96.25, 97.5, 100, ...
                                105, 110,130,200,400,500, 575,595, 610, 620, 625, 630, 635 640, 655, 680,710, 730, 735, 740, 745, ...
                                755, 770, 790, 825, 900, 1000 ];
                
                plot_IV(params_oo_BAPTA,scale_factors, 'g'); hold on;
                plot_IV(params_oo_EGTA,scale_factors, 'b');    
                plot_IV(params_wc_BAPTA,scale_factors, 'r');
                plot_IV(params_wc_EGTA,scale_factors, 'c');
                plot_IV(params_oo_NOBUFFER,scale_factors, 'k'); hold on;
         
                axis([0.1 600 0 0.04])
                xlabel('|I_{Ca++}| (pA)'); set(gca,'XTick',[0.1 1 10 100 1000],'XTickLabel',{'0.1','1','10','100', '10^3'});
                ylabel('\rho_{Ca++} (M)'); set(gca,'YTick',[1e-9,1e-8,1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]);

       % Save Figure             
       save_figure('Fig_Buffering_Calcium_V2',[10 10*4/3]);

       
       
       
function params = initialize_state(config, buffer)
        
	% Config Options - 'inside-out','outside-out','whole-cell'
    % Buffer Options - 'BAPTA', 'EGTA'
    
	% Solution Parameters
        % Ions are  Ca++,           EGTA/BAPTA,   Ca-EGTA/BAPTA,    K+,    Cl- 
    	params.D =  [0.4048,     0.25,    0.25,          1,    1.0382]' * 1.957e-9; % Diffusion Coefficient
                    % Ca++ - Barry and Lynch, J Membr Biol (1991) Vol 121, p101
                    % EGTA - Table 2 - Keramidas et. al. Journal of Neuroscience Methods (1999) Vol 89, p41–47. 
                    % BAPTA - Commonly assumed to match EGTA 
                    
        % Buffer Rate Constants 
            % Naraghi, Mohammad. 1997. “T-Jump Study of Calcium Binding Kinetics of Calcium Chelators.” 
            % Cell Calcium 22 (4): 255–68 -> pH 7.2
        Kon_EGTA  = 2.7e3; % Assocation rate of Ca++ and EGTA in mM/second
        Koff_EGTA = 0.5 ;  % Dissociation rate of Ca-EGTA in seconds
        Kon_BAPTA = 4.5e5  ;  % Assocation rate of Ca++ and BAPTA in mM/second
        Koff_BAPTA = 79 ; % Dissociation rate of Ca-BAPTA in seconds
        
    
        % Buffer parameters
        switch buffer
            case 'BAPTA'
                params.Kon =  Kon_BAPTA  ;  % Assocation rate 
                params.Koff = Koff_BAPTA ;  % Dissociation rate 
                params.Kd = Koff_BAPTA / Kon_BAPTA; % Equilibrium constant.      
                % Ions are  Ca++,  BAPTA,   Ca-BAPTA,  K+,    Cl- 
                params.z =  [2,     -4,    -2,         1,      -1]';  % Relative charge for BAPTA at pH 7
                
            case 'EGTA'    
                params.Kon =  Kon_EGTA  ;  % Assocation rate 
                params.Koff = Koff_EGTA ;  % Dissociation rate 
                params.Kd = Koff_EGTA / Kon_EGTA; % Equilibrium constant.      
                % Ions are  Ca++,  EGTA, Ca-EGTA,     K+,    Cl- 
                params.z =  [2,     -2,    0,         1,      -1]';  % Approximately right for EGTA at pH 7
                                                                     % See Keramides   
        end
        
        % Solutions
        params.rho_pip = [ 0 , 10, 0, 150, 150 ]';  % 150mM KCl in pipette + 10mM EGTA or BAPTA 
        idx = [1:4]; params.rho_pip(5) = sum(params.rho_pip(idx) .* params.z(idx));
        params.rho_bath = [2, 0, 0, 150, 154]';    % 150mM KCl and 2mM CaCl2  in bath
        
	% Pipette parameters
        params.r_tip = 0.31864e-6; 	% Tip radius of pipette (in m)
        params.theta =  10 * pi/180;
        params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
        params.r_pip = 0.5e-3 ; 	% Pipette inner radius0. of 0.5mm at the pipette electrode
		
        params.V_hold = -60e-3; 	% Initial holding voltage
        %params.Rseries = 0e6;       % Series resistance correction in Ohms
        %params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                    % Must be greater than Rpip * Cpip for
                                    % stability.
                                    
	% Membrane parameters
        params.membrane_func = @membrane_current_GHK; 	% Function describing membrane permeability to ions
		params.membrane_Pf = 0e-5   ; 			% Membrane permeability to water (m/s)
		
    % Initialize segments
      switch config
            
            case 'outside-out'
        		params.Nseg  =  61; 		% Number of segments between the pipette electrode and tip
        		params.membrane_x = 1e-9;   			% Position of membrane patch back from the tip (in m)
                params.membrane_perm = [1.74e-06 0 0 0 0]' / 1.0232  ; 		% Membrane permability to each ion in m/s
        
                
                % Custom segment positions to make equispaced on log scale
                 xmin = -1 *  params.r_pip / tan( params.theta/2);                % Pipette position
                 xmax = -1 * params.membrane_x;                                  % Membrane position
                 xmax = -1 * params.r_tip / tan( params.theta/2) - params.membrane_x;
                 params.membrane_N = params.Nseg - 1;                            % Index of membrane
                 idx = [1:params.Nseg];                                          % Index of positions
                 mu = 0.91; % Ratio of distance between successive segmetns (e.g. spacing of 1, 0.8, 064, ...)
                 pos = 1-cumsum(mu.^(idx-1)); pos = pos/pos(params.membrane_N);% Scaled position of each segment
                 params.x = 1 ./ ( 1/xmin + (1/xmax - 1/xmin) * pos) - xmax - params.membrane_x;
                 params.x(params.Nseg) = 0 ;% Segment positions
                 %fprintf('%f\n',params.x*1e9);
                 

            case 'inside-out'
                params.config = 'inside-out';
                rho_pip =  params.rho_pip ;
                rho_bath = params.rho_bath;
                params.rho_pip  = rho_bath;
                params.rho_bath = rho_pip;
                
                params.membrane_x = 4e-6;   % Membrane positioned 6 microns mm back from tip
                N_intra = 40; % Number of segments
                N_extra = 10; % Regulat number of segments
                
                % Logarithmic spacing for intra-cellular side
                mu = 0.8; % Ratio of distance between successive segmetns (e.g. spacing of 1, 0.8, 064, ...)
                idx = [1:(N_intra-1)];
                sf = [0, cumsum(mu.^(idx-1))];
                xi = params.membrane_x * sf / max(sf);
                
                % 1/r for extra-cellular
                um = 1/ params.membrane_x;
                ua = tan( params.theta/2) / params.r_pip;  
                
                xo = 1./ ([0:(N_extra-1)]/(N_extra-1) * (ua-um) + um );
            
                params.x = [xi, xo]; % Position of segmemts
                params.membrane_N = N_intra; % Position of membrane
                params.Nseg  =  length(params.x); 		% Number of segments between the pipette electrode and tip
                params.membrane_perm = [1.74e-06 0 0 0 0]' /6.9675 / 1.0232  ; 		% Membrane permability to each ion in m/s
                
            case 'whole-cell'

                % Cell and Membrane parameters
            	params.A_cell = 1e-9; % Cell Area
                params.Lambda_Cell = 1e-15 ; % Cell volume
            
            % Calculate positions and segments
                Npip  = 11 ; % Number of pipette segments
                Ncell = 5  ; % Number of segments from pipette to cell surface
                Nmem  = 12 ; % Number of segments at membrane
                mem_width = sqrt( params.D(1) / (params.Kon * params.rho_pip(2)) ) * 5;
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
                router = Rcell - mem_width;
                xcell = [ 1./ (1/rinner + (1/router-1/rinner)*[1:Ncell]/Ncell), ... %  1/r spacing for cell
                      router + mem_width * [1:Nmem]/Nmem]-rinner;                    %  mem_width/Nmem for region close to membrane
                Acell = omega_cell * (xcell+rinner).^2;
            
            % Bath segments
                rmem = Rcell + mem_width;
                rbath = params.r_pip;
                xbath = [ Rcell + mem_width *[0:Nmem]/Nmem, ...     % mem_width/Nmem for region close to membrane
                          1./ (1/rmem + (1/rbath - 1/rmem)* [1:Ncell]/Ncell) ] - rinner; % 1/r spacing out to bath pipette
                Abath = omega_cell * (xbath+rinner).^2;
    
            % Combine segments.
            params.x = [xpip, xcell, xbath];
            params.A = [Apip, Acell, Abath];
            params.membrane_N = Npip + Ncell + Nmem; % Position of membrane
            params.Nseg  =  length(params.x); 		% Number of segments between the pipette electrode and tip
            params.membrane_perm = [1.74e-06 0 0 0 0]' /3240 / 1.0232 ; 		% Membrane permability to each ion in m/s
      
      end
        
        
        % Initialize all remaining parameters and buffering.
        params.Constraints.type = 'buffer';
    
        % OO
        params = initialize_system(params);
        params_oo_BAPTA.Constraints.err_tol{5}  = 1e-7; % rate error in mM /second    
        params_oo_BAPTA.Constraints.stepsize{5}  = 1e-4; % step size  in mM /second    
       

 function plot_rho_vs_x_ss(params,plot_options)
        
        % Find SS Conditios
        [params, results ] = evolve_backward_euler(params, inf );
        [paramst, results ] = evolve_backward_euler(params, inf );

        fprintf(' %f pA',results.I_pip*1e12); % report Current
        
        
        % Plot Current
        subplot('Position', plot_options.panels(1,:));
            idx = [1:(paramst.Nseg-1)]; x = paramst.x(idx) * -1e6;
            semilogx(x, abs(paramst.J(1,:) * paramst.F * paramst.z(1))*1e12,'k'); hold on;
            semilogx(x, abs(paramst.J(3,:) * paramst.F * paramst.z(1))*1e12,'g'); hold on;
            axis(plot_options.axes(1,:));
            h = gca; set(gca,'TickDir','out','XTick',[]);
            ylabel('|I| (pA)');
            
        subplot('Position', plot_options.panels(2,:));
            x = paramst.x * -1e6;
            loglog(x, paramst.rho(1,:)*1e-3,'k'); hold on;
            loglog(x, paramst.rho(2,:)*1e-3,'b'); hold on;
            loglog(x, paramst.rho(3,:)*1e-3,'g'); hold on;
            axis(plot_options.axes(2,:));
            h=gca;
            %set(h,'TickDir','out'); 
            set(h,'XTick',[]);
            ylabel('\rho (M)');set(h,'YTick',10.^[-9:2:-2]);
            
        subplot('Position', plot_options.panels(3,:));
            loglog(x, paramst.rr(3,:)/1e3, 'g'); hold on;
            axis(plot_options.axes(3,:));
            h=gca;
            ylabel('d\rho_{Ca++-buffer}/dt (M/s)'); 
            set(h,'YTick',10.^[-3:2:1]);
            xlabel('|x| (\mum)'); set(h,'XTick',10.^[-3:1]);
            

     function plot_IV(params,scale_factors, plotcol)
         % Calculate and plot SS IV
         
            perm = params.membrane_perm;
         
            for j=1:length(scale_factors)
                params.membrane_perm = perm * scale_factors(j); % Membrane Permeability
                [params, results ] = evolve_backward_euler(params, inf );
                I(j) = abs(results.I_pip)*1e12;
                rho(j) = results.rho_intra(1) * 1e-3;
                %fprintf(1,' %f, ',scale_factors(j));
            end
   
            loglog( I, rho, plotcol);
            