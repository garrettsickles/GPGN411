% Plot Tensor Components

% First we must import the data
T_xx = csvread('T_xx_m.csv');
T_xy = csvread('T_xy_m.csv');
T_xz = csvread('T_xz_m.csv');
T_yy = csvread('T_yy_m.csv');
T_yz = csvread('T_yz_m.csv');
T_zz = csvread('T_zz_m.csv');

x_obs = csvread('x_obs.csv');
y_obs = csvread('y_obs.csv');


% Then we must plot it!
subplot(3, 3, 1);
surfc(x_obs, y_obs, T_xx, 'EdgeColor', 'none');
title('Txx');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'Eotvos');

subplot(3, 3, 5);
surfc(x_obs, y_obs, T_yy, 'EdgeColor', 'none');
title('Tyy');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'Eotvos');

subplot(3, 3, 9);
surfc(x_obs, y_obs, T_zz, 'EdgeColor', 'none');
title('Tzz');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'Eotvos');

subplot(3, 3, 2);
surfc(x_obs, y_obs, T_xy, 'EdgeColor', 'none');
title('Txy');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'Eotvos');

subplot(3, 3, 6);
surfc(x_obs, y_obs, T_yz, 'EdgeColor', 'none');
title('Tyz');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'Eotvos');

subplot(3, 3, 3);
surfc(x_obs, y_obs, T_xz, 'EdgeColor', 'none');
title('Txz');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'Eotvos');