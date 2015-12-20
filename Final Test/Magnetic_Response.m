function Magnetic_Response()
    clear;
    system('./finalTest');
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
    
    [I, D] = Helbig_Method(E_obs, N_obs, B_x, B_y, B_z);
    disp([num2str(I),', ',num2str(D)]);
end

function [Inc, Dec] = Helbig_Method(X, Y, Bx, By, Bz)
    fBx = @(x,y) x.*interp2(X,Y,Bz,x,y,'cubic');
    fBy = @(x,y) y.*interp2(X,Y,Bz,x,y,'cubic');
    fBz = @(x,y) x.*interp2(X,Y,Bx,x,y,'cubic');
    mx = (-1./(2*3.14159))*integral2(fBx,min(X),median(X),min(Y),median(Y),'RelTol',1e-4,'AbsTol',1e-6);
    mx = mx + (-1./(2*3.14159))*integral2(fBx,median(X),max(X),median(Y),max(Y),'RelTol',1e-4,'AbsTol',1e-6);
    my = (-1./(2*3.14159))*integral2(fBy,min(X),median(X),min(Y),median(Y),'RelTol',1e-4,'AbsTol',1e-6);
    my = my + (-1./(2*3.14159))*integral2(fBy,median(X),max(X),median(Y),max(Y),'RelTol',1e-4,'AbsTol',1e-6);
    mz = (-1./(2*3.14159))*integral2(fBz,min(X),median(X),min(Y),median(Y),'RelTol',1e-4,'AbsTol',1e-6);
    mz = mz + (-1./(2*3.14159))*integral2(fBz,median(X),max(X),median(Y),max(Y),'RelTol',1e-4,'AbsTol',1e-6);
    
    m =(mx.^2 + my.^2 + mz.^2).^(1/2);
    Inc = asind(mz./m);
    Dec = atand(my./mx);
end