% Plot Magnetic Field Components

% First we must import the data
B_x = csvread('B_x_m.csv');
B_y = csvread('B_y_m.csv');
B_z = csvread('B_z_m.csv');
Total = csvread('Total_m.csv');

N_obs = csvread('x_obs.csv');
E_obs = csvread('y_obs.csv');

% Then we must plot it!
subplot(2, 2, 1);
surfc(E_obs, N_obs, B_x, 'EdgeColor', 'none');
view(2);
title('B_x');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'nT');

subplot(2, 2, 2);
surfc(E_obs, N_obs, B_y, 'EdgeColor', 'none');
view(2);
title('B_y');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'nT');

subplot(2, 2, 4);
surfc(E_obs, N_obs, B_z, 'EdgeColor', 'none');
view(2);
title('B_z');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'nT');

subplot(2, 2, 3);
surfc(E_obs, N_obs, Total, 'EdgeColor', 'none');
view(2);
title('Total');
xlabel('Easting (m)');
ylabel('Northing (m)');
colorbar('location','eastoutside');
c = colorbar('location','eastoutside');
xlabel(c, 'nT');