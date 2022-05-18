% Demonstrate that Equivalent circuit can lead to logical inconsistencies
%
% Scenario is a whole-cell containing 150mM KCl held at 100mV
%
% At t = 0 membrane becomes permeable to K+ with gK+ = 20nS -> 2nA K+
% outward current
%
% Plot I, V, drho_dt for Equivalent Electric circuit
% 

   g_mem_K = 20e-9;        % Membrane conductance in Siemens - 2nA initial current
   Cmem = 10 * 1e-12;       % Membrane conductance in pF
   Lambda_cell = 1e-15;     % Cell volume in m^3
   F = 96485 ; % Faraday's constant
   
   Vpip = 100 * 1e-3 ;      % Holding voltage of pipette in V
   Rpip = 5 *1e6;           % Membrane resistance in Ohm
   rho = 150;               % Initial concentration of KCl in mM
   D_K  = 1.957e-9 ;            % K+ diffusion coefficient
   D_Cl = 1.0382* D_K;   % Cl-1 diffusion coefficient
   g_pip_K =  D_K  / (D_K + D_Cl) * 1/Rpip; % Pipette conductance
   g_pip_Cl = D_Cl / (D_K + D_Cl) * 1/Rpip; % Pipette conductance
   
   %Timebase
   t_val_q = 1e-6 * [-20:230]; % Timescale for electric steady-state
   t_val_rho = 1e-5* [24:17e5]; % Timescale for concentration changes
   t_val = [t_val_q, t_val_rho];
   Nq = length(t_val_q);
   idx_q = [1:Nq];
   Ne = length(t_val);
   idx_rho = [(Nq+1):Ne];

   I_mem = zeros(size(t_val)); % Membrane Current
   g_mem = (t_val >= 0) * g_mem_K; % Membrane conductance
   V_cell = Vpip * ones(size(t_val));% Cell voltage
   I_pip_K = zeros(size(t_val));% Pipette current
   I_pip_Cl = zeros(size(t_val));% Pipette current
   drho_K_dt = zeros(size(t_val));% Rate of change of [K+]cell
   drho_Cl_dt = zeros(size(t_val));% Rate of change of [Cl-]cell
   rho_K = rho * ones(size(t_val));% [K+]cell
   rho_Cl = rho * ones(size(t_val));% [Cl-]cell
   
   
   
   % Calculate evoluation using simple forward Euler method.
   for j=2:length(t_val)
       I_mem(j) = V_cell(j-1) * g_mem(j);
       I_pip_K(j) = (Vpip - V_cell(j-1)) * g_pip_K;
       I_pip_Cl(j) = (Vpip - V_cell(j-1)) * g_pip_Cl;
       V_cell(j) = V_cell(j-1) + (I_pip_K(j-1) + I_pip_Cl(j-1) - I_mem(j-1)) * ( t_val(j) - t_val(j-1) ) / Cmem;
       drho_K_dt(j) = (I_pip_K(j) - I_mem(j))/ (F * Lambda_cell);
       drho_Cl_dt(j) = (-I_pip_Cl(j))/ (F * Lambda_cell);
       rho_K(j) = rho_K(j-1) + drho_K_dt(j) * ( t_val(j) - t_val(j-1) );
       rho_Cl(j) = rho_Cl(j-1) + drho_Cl_dt(j) * ( t_val(j) - t_val(j-1) );
   end
   I_elec = ( I_pip_K + I_pip_Cl );
   
   
   
   %Initialize Figure
   figure(1); clf;
   fig_dim = [12 10];
   fig_name = 'Fig_Equiv_Circs_illogical_v2';
   save_figure(fig_name,fig_dim);
   
   
   % Position of rows and columns for figure.
        yp = [0.08, 0.37, 0.57, 0.75]; % Position of rows
        yh = [ 0.25, 0.15,  0.15, 0.2]; % Height of rows
        xp = [0.1 0.5]; % Position of columns
        xw = [0.375 0.375];
        t_range_q = [-20 230];
        t_range_rho = [0 17];
   
        
        
   % Plot Currents on Top
   subplot('position', [xp(1) yp(1) xw(1), yh(1) ]);
   y_range = [0 2];
   tbase = t_val*1e6;
   
   
   % Plot Currents on Top
   subplot('position', [xp(1) yp(4) xw(1), yh(4)]);
   y_range = [0 2];
   tbase = t_val_q*1e6;
    plot(tbase, I_elec(idx_q)*1e9 ,'k--'); hold on; 
    plot(tbase, I_mem(idx_q)*1e9,'k'); 
    plot(tbase, 1e9 * I_pip_K(idx_q),'b'); % K+ Pipette Current 
    plot(tbase, 1e9 * I_pip_Cl(idx_q),'g'); % Cl- Pipette Current
   axis([t_range_q, y_range]);
   box off; h = gca;  h.XAxis.Visible = 'off';    set(h,'TickDir','out');
   ylabel('I (nA)'); 
        
   subplot('position', [xp(2) yp(4) xw(2), yh(4)]);
   tbase = t_val_rho;
   plot(tbase, 1e9 * I_mem(idx_rho) ,'k'); hold on; 
   plot(tbase, 1e9 * I_elec(idx_rho),'k--'); 
   plot(tbase, 1e9 * I_pip_K(idx_rho),'b');
   plot(tbase, 1e9 * I_pip_Cl(idx_rho),'g');
   box off ; h = gca; h.YAxis.Visible = 'off';   h.XAxis.Visible = 'off';    
   axis( [t_range_rho, y_range]);
 
   % Plot Voltage Next
   idx_v = 3;
   subplot('position', [xp(1) yp(idx_v) xw(1), yh(idx_v)]);
   y_range = [90 100];
        
   tbase = t_val_q*1e6;
   plot(tbase,V_cell(idx_q)*1e3,'k'); 
   axis([t_range_q, y_range ]);
   box off; h = gca;  h.XAxis.Visible = 'off';    
   ylabel('V_{cell} (mV)');set(h,'TickDir','out');
        
   subplot('position', [xp(2) yp(idx_v) xw(2), yh(idx_v)]);
   tbase = t_val_rho;
   plot(tbase,V_cell(idx_rho)*1e3,'k');
   axis([t_range_rho, y_range ]);
   box off ; h = gca; h.YAxis.Visible = 'off'; h.XAxis.Visible = 'off';        
   
   % Plot drho_dt
   idx_v = 2;
   subplot('position', [xp(1) yp(idx_v) xw(1), yh(idx_v)]);
   y_range = [-22 0]; 
        
   tbase = t_val_q*1e6;
   plot(tbase, drho_K_dt(idx_q), 'b'); hold on;
   plot(tbase, drho_Cl_dt(idx_q),'g--'); 
   axis([t_range_q, y_range ]);
   box off; h = gca;  h.XAxis.Visible = 'off';  set(h,'TickDir','out');
   ylabel('d\rho_{cell}/dt (mM/s)');
   
   
   subplot('position', [xp(2) yp(idx_v) xw(2), yh(idx_v)]);
   tbase = t_val_rho;
   plot(tbase, drho_K_dt(idx_rho), 'b'); hold on;
   plot(tbase,drho_Cl_dt(idx_rho), 'g--'); 
   axis([t_range_rho, y_range ]);
   box off; h = gca; h.YAxis.Visible = 'off'; h.XAxis.Visible = 'off';   
      
   % Plot Concentration
   idx_v = 1;
   subplot('position', [xp(1) yp(idx_v) xw(1), yh(idx_v)]);
   y_range = [0 150.1]; 
        
   tbase = t_val_q*1e6;
   plot(tbase,rho_K(idx_q),'b'); hold on;
   plot(tbase,rho_Cl(idx_q),'g--'); 
   axis([t_range_q, y_range ]);
   box off; xlabel('(\mus)');
   ylabel('\rho_{cell} (mM)'); set(gca,'YTick',[0 50 100 150],'TickDir','out');
        
   subplot('position', [xp(2) yp(idx_v) xw(2), yh(idx_v)]);
   tbase = t_val_rho;
   plot(tbase,rho_K(idx_rho),'b'); hold on;
   plot(tbase,rho_Cl(idx_rho),'g--'); 
   axis([t_range_rho, y_range ]);
   box off; h = gca; h.YAxis.Visible = 'off';      
   xlabel('(s)'); set(gca,'XTick',[5 10 15],'TickDir','out'); 
   

   % Save the figure     
   save_figure(fig_name,fig_dim);
   
       
        
             
        
         
        
  