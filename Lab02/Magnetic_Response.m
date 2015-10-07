function Magnetic_Response(strength, inc, dec, kappa)
pi = 3.1415926;

% First we must import the data
T_xx = csvread('T_xx_m.csv');
T_xy = csvread('T_xy_m.csv');
T_xz = csvread('T_xz_m.csv');
T_yx = (-1.0) .* T_xy;
T_yy = csvread('T_yy_m.csv');
T_yz = csvread('T_yz_m.csv');
T_zx = (-1.0) .* T_xz;
T_zy = (-1.0) .* T_yz;
T_zz = csvread('T_zz_m.csv');

northing = csvread('x_obs.csv');
easting = csvread('y_obs.csv');

B_n(1,1) = strength * cos(inc * pi / 180) * cos(dec * pi / 180);
B_n(2,1) = strength * cos(inc * pi / 180) * sin(dec * pi / 180);
B_n(3,1) = strength * sin(inc * pi / 180);
B_nhat = transpose(B_n) ./ strength;

[m,n] = size(T_xx);

%Establish the total field anomaly
Total = zeros(m, n);
B_x = zeros(m, n);
B_y = zeros(m, n);
B_z = zeros(m, n);

%Iterate over all of the observation points and calculate the total-field
%   anomaly at each point
for i=1:m
    for j=1:n
        T = (kappa/(4*pi)) .* ([T_xx(i,j) T_xy(i,j) T_xz(i,j); T_yx(i,j) T_yy(i,j) T_yz(i,j); T_zx(i,j) T_zy(i,j) T_zz(i,j);] * B_n);
        B_x(i,j) = T(1);
        B_y(i,j) = T(2);
        B_z(i,j) = T(3);
        Total(i,j) = (B_nhat * T);
    end
end

%Plot the total-field anomaly as a interpolated 2D plot
figure;
surf(northing, easting, Total,'EdgeColor', 'None', 'facecolor', 'interp');
view(2);
title('Total Field Anomaly');
xlabel({'Easting'});
ylabel({'Northing'});
cb = colorbar('location','eastoutside');
xlabel(cb, 'nT');

figure;
surf(easting, northing, B_x,'EdgeColor', 'None', 'facecolor', 'interp');
view(2);
title('B_x');
xlabel({'Easting'});
ylabel({'Northing'});
cb = colorbar('location','eastoutside');
xlabel(cb, 'nT');

figure;
surf(easting, northing, B_y,'EdgeColor', 'None', 'facecolor', 'interp');
view(2);
title('B_y');
xlabel({'Easting'});
ylabel({'Northing'});
cb = colorbar('location','eastoutside');
xlabel(cb, 'nT');

figure;
surf(easting, northing, B_z,'EdgeColor', 'None', 'facecolor', 'interp');
view(2);
title('B_z');
xlabel({'Easting'});
ylabel({'Northing'});
cb = colorbar('location','eastoutside');
xlabel(cb, 'nT');
clear;
end