%
function Fig_Hv1_Giant_Patch_vs_WholeCell


    % I/pH versus T for voltage step protocol - pH 5.5, pH 7.5
    figure(1); clf;
    
    simulate_voltage_step('pH=5.5',0);
    simulate_voltage_step('pH=7.5',1);
    save_figure('Fig_Hv1_Macropatch_vs_wholecell_I_pH_v1',[ 10 12]);  
    
    
    % Detect pH change with ramp
    figure(2);clf
    simulate_voltage_step_with_ramp('pH=7.5');
    save_figure('Fig_Hv1_Macropatch_vs_wholecell_rampt',[ 10 6]);  
                
   % Plot I versus v - Column 1 for pH 5.5, Column 2 for pH 7.5
   % figure(3); clf;
   % calculate_iv('pH=5.5',1)
   % calculate_iv('pH=7.5',2)


 function simulate_voltage_step(solution, offset)

    V_hold = -80e-3; % Holding Voltage
    V_stim = 140e-3; % Stimulus Voltage
    V_tail = -80e-3; % Voltage for measuring tail current
    intervals.duration =            [1,         0.5,        4.5 ,    0.1,       0.3,    0.6   ];    % Time intervals in protocol (seconds)
    intervals.N_eval =              [2,         10,          20,      10,         10,       5];      % Number of time points during each time interval
    intervals.variables.V_hold =    {V_hold,    V_stim,     V_stim,  V_tail,    V_tail, V_tail};   % Voltage during each interval in V
    protocol = create_protocol(intervals);
    results.t = protocol.t;
    
    % Three configs we will simulate
    config = {'ideal', 'inside-out', 'whole-cell'};
    
    % Panel Positions
    l1 =0.05; w = 0.26; l2 = 0.38 ; l3 = 0.71;
    b1 = 0.03; h1 = 0.2; h2 = 0.2; 
    pos      = [l1, 0.5+2*b1+h1, w, h2; % pH5.5 control - I 
                l1, 0.5+b1,      w, h1; % pH 5.5 control- pH
                l2, 0.5+2*b1+h1, w, h2; % pH 5.5 patch - I 
                l2, 0.5+b1,      w, h1; % pH 5.5 patch- pH  
                l3, 0.5+2*b1+h1, w, h2; % pH 5.5 whole cell - I 
                l3, 0.5+b1,      w, h1; % pH 5.5 whole cell - pH  
                l1, 2*b1+h1, w, h2;     % pH 7.5 control - I 
                l1, b1,      w, h1;     % pH 7.5 control- pH
                l2, 2*b1+h1, w, h2;     % pH 7.5 patch - I 
                l2, b1,      w, h1;     % pH 7.5 patch- pH  
                l3, 2*b1+h1, w, h2;     % pH 7.5 whole cell - I 
                l3, b1,      w, h1;];     % pH 7.5 whole cell - pH  

 for j=1:length(config)

        params = init_config(solution,config{j},'dynamic');
        params.V_hold = V_hold; [params, res ] = evolve_backward_euler(params, inf); % Initialize State
        [params_f, res] = evolve_backward_euler(params, protocol);
        I = res.I_pip * 1e12; % Current in pA
        pH_intra = 3 - log10(res.rho_intra(1,:));
        pH_extra = 3 - log10(res.rho_extra(1,:));
    
        % I versus t
        subplot('position',pos(offset*6 +j*2-1,:));
            plot(results.t, I, 'k'); hold on
                axis([0 7 -410 700]);    axis off
                if (j==1)
                    plot(1 +  [0 2 ],-410 + [0 0],'k','Linewidth',2); % 2s scale bar
                    plot(1 +  [0 0] , 100 + [0 500],'k','Linewidth',2); % 500pA scale bar
                end
                
        % pH versus t    
        subplot('position',pos(offset*6 +j*2,:)); 
            plot(results.t, pH_intra, 'k'); hold on
            plot(results.t, pH_extra, 'k:');
            if (offset==0)

                 axis([0 7 5.39 5.8]); 
                 box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out';
                 ylabel('pH');set(gca,'YTick',[5.4:0.1:5.7]);
            else     
                axis([0 7 7.39 7.8]); 
                box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out';
                ylabel('pH');set(gca,'YTick',[7.4:0.1:7.7]);
                 
            end
            
    end
    
function simulate_voltage_step_with_ramp(solution)

    V_hold = -80e-3; % Holding Voltage
    V_stim = 140e-3; % Stimulus Voltage
    V_tail = -80e-3; % Voltage for measuring tail current
    V_ramp_up = [-80e-3, 140e-3]; % Ramp up
    V_ramp_d= [140e-3, -80e-3]; % Ramp down
    t_ramp = 0.22; % 1 V/s ramp speed.
    N_ramp = 21;
    
    % Short Protocol
    intervals.duration =            [0.78,         t_ramp,     t_ramp,     0.78,      0.5,        1.5 ,   t_ramp,       0.2,    0.58   ];    % Time intervals in protocol (seconds)
    intervals.N_eval =              [2,         N_ramp,     N_ramp,     2,          5,          6,   N_ramp,         5,       5];      % Number of time points during each time interval
    intervals.variables.V_hold =    {V_hold,    V_ramp_up,  V_ramp_d,   V_hold,V_stim,     V_stim,  V_ramp_d,    V_tail, V_tail};   % Voltage during each interval in V
    protocol_short = create_protocol(intervals);
    
    % Long Protocol
    intervals.duration =            [0.78,         t_ramp,     t_ramp,     0.78,      0.5,        7.5 ,   t_ramp,       0.2,    0.58   ];    % Time intervals in protocol (seconds)
    intervals.N_eval =              [2,         N_ramp,     N_ramp,     2,          5,         30,   N_ramp,         5,       5];      % Number of time points during each time interval
    intervals.variables.V_hold =    {V_hold,    V_ramp_up,  V_ramp_d,   V_hold,V_stim,     V_stim,  V_ramp_d,    V_tail, V_tail};   % Voltage during each interval in V
    protocol_long = create_protocol(intervals);
    idx_r1 = 3;
    idx_r2 = 7;

    protocols = {protocol_short, protocol_long};
    t_range = { [0 5], [0 11]};
  
    % Plot positions
    le1= 0.1; wid1 = 0.1615; % column 1 Short stim
    le2= 0.3115 ; wid2 =  0.3551  ; % Column 2 Long stim
    le3 = 0.67; wid2 = 0.31; % Colum 3 IV
    pos = [le1 0.85  wid1 0.12; % voltage protocol
           le1 0.41  wid1 0.4; % Current
           le1 0.1   wid1 0.3 ; % pH
           le2 0.85  wid2 0.12; % voltage protocol
           le2 0.41  wid2 0.4; % Current
           le2 0.1   wid2 0.3 ; % pH
           le3 0.2  wid2 0.7 ;]; % IV
       
    % Initialize State;
    params = init_config(solution,'whole-cell','dynamic');
    params.V_hold = V_hold;  [params, res ] = evolve_backward_euler(params, inf); % Initialize State
    
    for j=1:length(protocols)
    
        [params_f, res] = evolve_backward_euler(params, protocols{j});
        
        I = res.I_pip * 1e12; % Current in pA
        pH_intra = 3 - log10(res.rho_intra(1,:));
        pH_extra = 3 - log10(res.rho_extra(1,:));

        % Plot
            % V_command
            subplot('position', pos(3*j-2,:));
                plot(protocols{j}.t, res.V_command*1e3,'k');
                axis([t_range{j} -100 200]); axis off
                axis off;
                
            % Current
            subplot('position',pos(3*j-1,:));
                plot(protocols{j}.t, I, 'k');  hold on
                plot(1 +  [0 2 ],-100 + [0 0],'k','Linewidth',2); % 2s scale bar
                plot(1 +  [0 0] , 100 + [0 250],'k','Linewidth',2); % 250pA scale bar
                axis([t_range{j} -200 400]);
                axis off

            % pH
            subplot('position',pos(3*j,:));
                plot(protocols{j}.t, pH_intra,'k'); hold on
                plot(protocols{j}.t, pH_extra,'k'); 
                axis([t_range{j} 7.39 7.8]); 
                box off; h = gca;  h.XAxis.Visible = 'off';    h.TickDir='out';
                ylabel('pH');set(gca,'YTick',[7.4:0.1:7.7]);

         % Ramps       
         idx1 = [protocols{j}.index(idx_r1,1):protocols{j}.index(idx_r1,2)];
         idx2 = [protocols{j}.index(idx_r2,1):protocols{j}.index(idx_r2,2)];
         V_ramp_init{j} = res.V_command(idx1)*1e3; I_init{j} = I(idx1);
         V_ramp_end{j}  = res.V_command(idx1)*1e3; I_end{j}  = I(idx2);
    
    end
    
    % IV Plot
    subplot('position',pos(7,:));
            plot(V_ramp_init{1}, I_init{1},'k--'); hold on
            plot(V_ramp_end{1},I_end{1},'b'); hold on
            plot(V_ramp_end{2},I_end{2},'k'); hold on
            
            % Add on true Vrev at beginning of ramp.
            %idx3 = protocol.index(idx_r2,1) - 1;
            %Vrev = -1 * params.R * params.T / params.F * log(res.rho_intra(1,idx3)/res.rho_extra(1,idx3))
            %plot( Vrev * 1e3+ [0 0], [0 -40],'k');
            
            ax = gca; box off; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';set(gca,'TickDir','out');
            axis([-30 30 -70 70]);
            set(gca,'XTick',[-20 -10 10 20]);
            set(gca,'YTick',[-50 -25 25 50]);
            xlabel('V (mV)');
            ylabel('I (pA)');
            

function params = init_config(solutions, config, membrane)        
    
    % Solutions
    switch solutions
        case 'pH=7.5'
                        % H+,   B,     HB,     TMA+,   methane-sulfate 
                        % Buffer is HEPES
            params.D =  [4.759,   0.3,    0.3,    0.611,      0.664               ]' * 1.957e-9; % Diffusion Coefficient
            params.z =  [1,     -1,     0,      1,     -1]';  % Relative charge
            
            pKa = 7.5; % pKa of buffer
            params.Kon = 1e10 / 1e3 ;               % Assocation rate of H+ and B in mM^-1 * sec^-1
            params.Koff = 10^(3-pKa) * params.Kon ; % Dissociation rate of HB in s^-1
            params.Kd =   10^(3-pKa) ;                % Kd of HB buffer in mM -> pKa = 3-log10(params.K_HB) 
                
            params.rho_pip =  init_buffer(params, 7.5, 100, 60); %100mM buffer at pH 7.5, 60mM TMA-MeSO3
            params.rho_bath = init_buffer(params, 7.5, 100, 60); %100mM buffer at pH 7.5, 60mM TMA-MeSO3
            Hv1_perm = 1.2654e-10/2; %  Absolute membrane permeability to give 5nS conductance at 140mV
            
        case 'pH=5.5'
                        % H+,   B,     HB,     TMA+,   methane-sulfate 
                        % Bufer is MES
            params.D =  [4.7,   0.37,    0.37,    0.611,      0.664]' * 1.957e-9; % Diffusion Coefficient
            params.z =  [1,     -1,     0,      1,     -1]';  % Relative charge
            
            pKa = 6.15; % pKa of buffer
            params.Kon = 1e10 / 1e3 ;               % Assocation rate of H+ and B in mM^-1 * sec^-1
            params.Koff = 10^(3-pKa) * params.Kon ; % Dissociation rate of HB in s^-1
            params.Kd =   10^(3-pKa) ;                % Kd of HB buffer in mM -> pKa = 3-log10(params.K_HB) 
                
            params.rho_pip =  init_buffer(params, 5.5, 100, 60); %100mM buffer at pH 7.5, 20mM TMA-MeSO3
            params.rho_bath = init_buffer(params, 5.5, 100, 60); %100mM buffer at pH 7.6, 60mM TMA-MeSO3
            Hv1_perm = 9.28e-13/2; %  Absolute membrane permeability to give 5nS conductance at 140mV
            
    end
    
    % Patch-Clamp Config
    switch config
        case 'whole-cell'
            params.config = 'whole-cell';      % Options - 'inside-out','outside-out','whole-cell'
     
        % Pipette parameters
            params.r_tip = 0.812e-6; 	% Tip radius of pipette (in m)
            params.theta = 10 * pi/180;
            params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
            params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
            %params.Rseries = 0e6;       % Series resistance correction in Ohms
            %params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                        % Must be greater than Rpip * Cpip for
                                        % stability.
                                    
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

        case 'inside-out'
            params.config = 'inside-cell';      % Options - 'inside-out','outside-out','whole-cell'
     
        % Pipette parameters
            params.r_tip = 5e-6; 	% Tip radius of pipette (in m)
            params.theta = 36 * pi/180;
            params.d_ratio= 0.75;        % Inner PIpette diameter/outer pipette diameter
            params.r_pip = 0.5e-3 ; 	% Pipette inner radius of 0.5mm at the pipette electrode
            %params.Rseries = 0e6;       % Series resistance correction in Ohms
            %params.Rseries_tau = 2e-6;  % Series resistance compensation time
                                        % Must be greater than Rpip * Cpip for
                                        % stability.
                                    
    	% Cell and Membrane parameters
            params.membrane_x = 15e-6; % Position of membrane
        
        % Calculate positions and segments
            Nintra = 11; % Number of segments from bath to near membrane region
            Nmem  = 12 ; % Number of segments at membrane
            mem_width = sqrt( params.D(1) / (params.Kon * params.rho_pip(2)) ) * 5;
            Nextra = 11; % Number of segments from near membrane to pipette electode
            
            % Pipette segments
            omega = 4 * pi * sin(params.theta/4)^2; %   Solid angle
            dtip =  params.r_tip / tan(params.theta/2); % Distance from source to tip.
            dintra = dtip + params.membrane_x - mem_width; % Distance from source to intra-cellular side of membrane-region
            dextra = dtip + params.membrane_x + mem_width; % Deistrance from source to extra-cellular side of membrane region
            dpip =  params.r_pip / tan(params.theta/2); % Distance from source to pipette electrode.
            
            
            xintra = [ 1./ (1/dtip + (1/dintra - 1/dtip) *([0:(Nintra-1)]/Nintra )) - dtip, ... % 1/r spacing from bath to near membrane
                       params.membrane_x + mem_width * [-Nmem:0]/Nmem; ];                      % Linear spacing near membrane
            xextra = [ params.membrane_x + mem_width * [0:Nmem]/Nmem, ...   % Linear spacing away from membrane
                       1./ (1/dextra + (1/dpip - 1/dextra) * [1:Nextra]/Nextra) - dtip ];  %1/r spacing to pipette electode
            params.x = [xintra, xextra];
            params.A = omega * (dtip + params.x).^2;
            params.membrane_N = Nintra + Nmem +1; % Position of membrane
            params.Nseg  =  length(params.x); 		% Number of segments between the pipette electrode and tip
        
      case 'ideal'
            params.config = 'ideal';      % Options - 'inside-out','outside-out','whole-cell'
            params.Nseg = 2;              % Number of segments
            
        % Calculate positions and segments
            params.x = [-1, 1 ] * 1e-6;
            params.A = 1e-9 * ones(size(params.x));
            params.membrane_N = 1; % Position of membrane
            params.Nseg  =  length(params.x); 		% Number of segments between the pipette electrode and tip       
     end
            

    % Common Properties
            params.Constraints.type = 'buffer';     
            params.V_hold = 0e-3; 	% Initial holding voltage
    
    % Membrane Properties        
            params.membrane_func = @membrane_current_Hv1; 	% Function describing membrane permeability to ions
        	params.membrane_perm = [ 0 0 0 0 0]'; 		% Membrane permability to each ion in m/s independent of Hv1
            params.membrane_Pf = 0e-5   ; 			% Membrane permeability to water (m/s)
            
            % Hv1 Paramters from table S1 of Villalba-Galea 2014
            % State 0 = C1
            % State 1 = C2
            % State 2 = C3
            % State 3 = O
            % State 4 = D - Hypothetical deactivation intermediate
            % Activation pathway
            a01 = 3.7e-4;   b10 = 3.6e-6;   za01 = 0.16;    zb10 = 2.42;      % Transitions between states 0 and 1
            a12 = 3.9e-5;   b21 = 1.1e-1;   za12 = 0.92 ;   zb21 = 1.52;      % Transitions between states 1 and 2
            a23 = 6.6e-3;   b32 = 1.6e-3;   za23 = 3e-5;    zb32 = 2.7e-7;    % Transitions between states 2 and 3
            
            % Deactivation pathway effectively O-> D followed by very slow
            % D -> C1
            a34 = 1.3e-3;     b43 =  5e-3;    za34 = -0.65;    zb43 = -1.35 ; % Transitions betwen state 3 and hypothetical deactivation intermediate state 4
            a04 = 1e-10;    b40 = 10;        za04 =1.5;       zb40 = 1.5;     % Ultra-slow transitions between state 4 back to 0th stage for recovery.
        
            params.membrane_Hv1_rate =  [ 0  a01  0    0   a04;   b10  0  a12    0  0; 0   b21 0 a23  0;  0  0   b32  0   a34;   b40 0 0  b43   0 ]; % Transition rates
            params.membrane_Hv1_gamma = [ 0  0    0    0     0;  -0.5  0  0.7    0  0; 0  -0.7 0 0.2  0;  0  0  -0.1  0     0;     0 0 0  0     0 ];  % Coupling to pH  
                            % look at Figure 5 of paper
            params.membrane_Hv1_z =      [ 0 za01 0    0  za04; -zb10  0 za12    0  0; 0 -zb21 0 za23 0;  0  0  -zb32 0  za34; -zb40 0 0 -zb43 0 ]; % Charge for transition
            params.membrane_Hv1_rate = params.membrane_Hv1_rate * 1e3; % Convert to seconds        
            params = initialize_system(params);
            params.membrane_Hv1_permH = Hv1_perm /params.A(params.membrane_N+1); % Membrane proton permeability when Hv1 fully open.   

     switch membrane
         case 'static'
            params.membrane_Hv1_state = [0 0 0 1 0]; % Start Hv1 in Conducting State

         case 'dynamic'
            params.membrane_Hv1_state = [1 0 0 0 0]; % Start Hv1 in Closed State
        
            % Add in channel degrees of freedom
            if strcmp(config,'ideal')
                % No density or reaction rate evolution
               	params.Constraints.param =      {'J',                       'V',                'Jsol',             'membrane_Hv1_state'}; 
				params.Constraints.fhandle =  	{@constraint_J,             @constraint_V, 		@constraint_Jsol,   @constraint_Hv1_State};
    			params.Constraints.segments =   {[1,params.Nseg-1],         [1,params.Nseg],	[1,params.Nseg-1],  params.membrane_N};
    			params.Constraints.components = {[1,size(params.rho,1)],    [1 1],              [1 1],              [1, length(params.membrane_Hv1_state)]};
    			params.Constraints.stepsize = 	{0.01,                      1e-5,               1e-18,              0.01};
                params.Constraints.err_tol = 	{1e-18,                     1e-5,               1e-19,              1e-6};

            else
                nchan = 6;
                params.Constraints.param{nchan}     = 'membrane_Hv1_state'; 
                params.Constraints.fhandle{nchan}   = @constraint_Hv1_State;
                params.Constraints.segments{nchan}  = params.membrane_N;
                params.Constraints.components{nchan}= [1 length(params.membrane_Hv1_state)];
                params.Constraints.stepsize{nchan} =  0.01;
                params.Constraints.err_tol{nchan}  =  1e-6;
                params = initialize_system(params);
            end
            
            % Add in channel degrees of freedom
            params = initialize_system(params);
            
     end   
    
     % Remove capacitance transients.   
        params.Ct = 0 * params.Ct;
        params.Cl = 0 * params.Cl;
    
     % Report Access Resistance and Series Resistance Correction
     Raccess = sum(resistance(params));
        fprintf('\n Access Resistance = %.2f MOhms\n', Raccess/1e6);
     
function result = init_buffer(params, pH, C_buf, C_NaCl)

        C_H = 1e3  / 10^pH; % Proton concentration in mM
        C_B = params.Kd / (params.Kd + C_H) * C_buf; % Free buffer
        C_HB = C_buf - C_B; % Bound buffer
        C_Na = C_NaCl + C_B - C_H; % Na+ concenrtation in (mM)
        C_Cl = C_NaCl ; % Cl- concentration (mM)
    
    	result = [ C_H , C_B, C_HB, C_Na, C_Cl   ]';  % Concentration in pipette


function calculate_iv(solution, column)
    
            config = {'inside-out', 'whole-cell'};
            linecol = {'k',     'b'};
            nrow =3 ; % 3 rows in figure.
            ncol = 2; % 2 columns.
    
            Vstep = 20e-3;
            Imax = 400; % Current range (in pA)
            
            for j=1:length(config)

                m = 0;
                
                for sign = -1:2:1
                    params = init_config(solution,config{j},'static');  % Initialize configuration
                    % Step up  voltage until reach max current       
                    Iabs  = 0;                
  
                    while (Iabs < Imax)
                    
                      params.V_hold = params.V_hold + Vstep * sign; 
                      [params, result] = evolve_backward_euler(params, inf); 
                      m=m+1; V(m) = params.V_hold; I(m) = result.I_pip*1e12; ;
                      pH_int(m) = 3-log10(result.rho_intra(1)); pH_ext(m) = 3-log10(result.rho_extra(1));
                      Iabs = abs(I(m));
                    end
                end
                
                
                % Sort results by voltage
                [V, idx] = sort(V);
                I = I(idx); pH_int = pH_int(idx); pH_ext = pH_ext(idx);   
                
                subplot(nrow,ncol,4+column)
                plot(I, pH_int,linecol{j}); hold on
                plot(I, pH_ext,[linecol{j},':']); hold on
            end
    
        if (column==1) 
            axis([-400 400 5 6]); % pH 5.5
        else
            axis([-400 400 7 8]); % pH 7.5
        end
        
        xlabel( 'I (pA)');
        ylabel('pH');
            
    
        
    
        