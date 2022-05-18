% function save_figure(fname, width_height)
%
%   Save figure as fname.eps, fname.jpg
%   If specified, width_height(1) = width of figure in cm
%                 width_height2) =  height of figure in cm
%


function save_figure(fname, width_height)

if (nargin<1)
    fname = 'figure';
end

if (nargin < 2)
    width = 12;
    height = 9;
else
    width = width_height(1); % Width of Figure in Inches
    height = width_height(2); % Height of Figure in Inches
end

% Set up the paper
% First set up the paper to be the right size.
fig = gcf;
set(fig,'PaperUnits','centimeters');
set(fig,'Units','centimeters');
papersize=get(fig,'PaperSize');
left = (papersize(1)-width)/2;
bottom= (papersize(2)-height)/2;
set(fig,'PaperPositionMode','manual');
mfs=[left,bottom,width,height];
set(fig,'PaperPosition',mfs);
set(fig,'Position',[0.0,0.0,width,height]);

% Finally save the figure 
eval(sprintf('print -depsc -tiff -r150 %s', sprintf('%s.eps',fname)));
saveas(fig, sprintf('%s.pdf',fname),'pdf');
eval(sprintf('print -djpeg -r300 %s',sprintf('%s.jpg',fname)));
eval(sprintf('print -dsvg -r300 %s',sprintf('%s.svg',fname)));
%plot2svg(sprintf('%s.svg',fname),fig,'jpg');