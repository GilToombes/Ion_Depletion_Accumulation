% function plot_Tau_Cell_vs_Raccess_Lambda
%
%   Plot dependence of Tau_Cell on access resistance and cell volume for
%   150mM KCl.

function fig_WholeCell_tau

    % Cell Volume
    Nl = 100; Lmin = 1e-15/20; Lmax = 1e-15 * 20;
	Lambda = Lmin * 10.^([0:Nl]/Nl * log10(Lmax/Lmin)) ;
    
    % Access Resistance
    Nr = 100; Rmin = 1e6; Rmax = 10e6;
    Raccess = Rmin * 10.^( [0:Nr]'/Nr * log10(Rmax/Rmin) );
    
    % Solution conductance.
    D = [1 , 1.0382]' * 1.957e-9; % Diffusion coefficient in m/s^2 for K+ and Cl-
    rho = [150, 150]';          % Concentraiton of each ion.
    R = 8.314; % Gas Constant
    F = 96485; % Faraday's Constant
	T = 293 ;  % Temperature
    sigma = sum(D.*rho) * F^2 / (R*T) ;

    % Calculate Time Course
    tau = Raccess * Lambda * sigma / D(1);
    clf; imagesc(log10(tau)); axis image;
    title('\tau_{cell} (s)');
    
    % Label X-axis
    XTick_val = [0.05, 0.2 1 5 20];
    XTick_Pos = interp1(Lambda, [1:length(Lambda)],XTick_val*1e-15  );
    set(gca,'XTick',XTick_Pos,'XTickLabel',XTick_val);
    xlabel(' Cell Volume - \Lambda (pL)');
    
    % Label Y-axis
    set(gca,'YDir','Normal');
    YTick_val = [1 2 5 10];
    YTick_Pos = interp1(Raccess, [1:length(Raccess)],YTick_val*1e6  );
    set(gca,'YTick',YTick_Pos,'YTickLabel',YTick_val);
    ylabel(' R_{access} (M\Omega)');
    
   
    % Colorbar
    Nc =256;
    colormap(jet(256));
    u = colorbar('SouthOutside');
    CTick_Val = [0.1 1 10 100];
    set(u,'XTick',log10(CTick_Val), 'XTickLabel', CTick_Val,'TickDir','out');
    
    set(gca,'TickDir','out');
    
    save_figure('Fig_Taucell_vs_R_Lambda_v1',[ 8 8]);       
        