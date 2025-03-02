function [] = setDefaultFigProp()
%Text size scaling set in save_figure
%Set some global parameters for general figure properties

%Set global properties for plotting for the executed script
set(groot, 'DefaultFigureColor', 'w'); %Default background color for figures

set(groot, 'DefaultAxesFontName', 'Arial'); %Font
set(groot, 'DefaultAxesLineWidth', 1.8); 
set(groot, 'DefaultLegendBackgroundAlpha', 0.5);

setappdata(groot, 'DefaultSubtitleFontSize', 10); %just a temp small subtitle font (save_figure should overwrite)

%Default Line
set(groot, 'DefaultLineLineWidth', 4); 
set(groot, 'DefaultLineMarkerSize', 16); 

%Default stair properties (inherit most from line)
set(groot, 'DefaultStairLineWidth', 4); 
set(groot, 'DefaultStairLineStyle', '-'); 

%Default scatter properties
%set(groot, 'DefaultScatterMarkerSize', 50);
%set(groot, 'DefaultScatterFilled', '1');
set(groot, 'DefaultScatterMarker', 'o'); 
set(groot, 'DefaultScatterSizeData', 200); % Set default marker size
set(groot, 'DefaultScatterLineWidth', 4);


end