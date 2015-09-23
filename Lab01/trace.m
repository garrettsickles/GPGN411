% Plot Trace Residuals (The Laplacian)

% First we must import the data
T_xx = csvread('T_xx_m.csv');
T_yy = csvread('T_yy_m.csv');
T_zz = csvread('T_zz_m.csv');

x_obs = csvread('x_obs.csv');
y_obs = csvread('y_obs.csv');

% Now we plot!
surfc(x_obs, y_obs, T_xx+T_yy+T_zz);
title('Residual Trace (Txx+Tyy+Tzz)');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'Eotvos');