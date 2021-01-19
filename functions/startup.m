%% Setting the default renderer to "Painters" instead of "OpenGL"
set(0, 'DefaultFigureRenderer', 'painters');

%% Setting font and size
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

%% Defining basic discrete colors based on brewermap

color_scheme = brewermap(9,'Set1');
color_scheme(10,:) = 115*ones(1,3)/255;


red = color_scheme(1,:);
blue = color_scheme(2,:);
green = color_scheme(3,:);
purple = color_scheme(4,:);
orange = color_scheme(5,:);
yellow = color_scheme(6,:);
brown = color_scheme(7,:);
pink = color_scheme(8,:);
% gray = color_scheme(9,:);
gray = [0.2 0.2 0.2];
darkgray = 115*ones(1,3)/255;
white = [1 1 1];

lightblue = [66,146,198]/255;
lightgreen = [161,217,155]/255;

%% Defining continuous color_scheme for colormap

color_scheme_OrRd = brewermap(100,'OrRd');
color_scheme_YlOrRd = brewermap(100,'YlOrRd');
color_scheme_BuPu = brewermap(100,'BuPu');

a = brewermap(100,'BuPu').';
b = fliplr(a(:,1:length(a)));
c = b.';
color_scheme_BuPuFlipped = c;

color_scheme_BuGn = brewermap(100,'BuGn');

color_scheme_BuBkWh(1,:) = [66,146,198]/255;
color_scheme_BuBkWh(2,:) = [0 0 0]/255;
color_scheme_BuBkWh(3,:) = [255,255,255]/255;

color_scheme_spectral9 = brewermap(9,'Spectral');
color_scheme_spectral9_mod = color_scheme_spectral9; %light blue (last)
color_scheme_spectral9_mod(9,:) = [237,248,251]/255;
% color_scheme_spectral9_mod(1,:) = [203,24,29]/255;

%% 

color_scheme_npg_hex = ['#E64B35';'#4DBBD5';'#00A087';'#3C5488';'#F39B7F';'#8491B4';'#91D1C2';'#DC0000';'#7E6148';'#B09C85'];
color_scheme_npg = hex2rgb(color_scheme_npg_hex);





