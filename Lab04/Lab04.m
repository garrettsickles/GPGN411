function Lab04()
    filename = uigetfile('*.dat','Input 2D data file');
    file = fopen(filename, 'rt');
    ny = fscanf(file,'%d',1);
    nx = fscanf(file,'%d',1);
    for j=1:ny,
        for i=1:nx,
            y(i,j) = fscanf(file,'%f',1);
            x(i,j) = fscanf(file,'%f',1);
            z(i,j) = fscanf(file,'%E',1);
        end
    end
    fclose(file);
    
    [Ky,Kx] = meshgrid((2*pi).*Fk(y(1,:)),(2*pi).*Fk(x(:,1)));
    
    Z = fftshift(fft2(z));
    amplitude = abs(Z);
    
    figure('position', [0, 0, 600, 450]);
    subplot(2,2,1);
    pcolor(y, x, z);
    xlabel('x (m)');
    ylabel('y (m)');
    title('Data');
    colorbar;
    colormap jet;
    shading interp;
    axis tight;

    subplot(2,2,2);
    pcolor(Ky, Kx, log10(amplitude.^2));
    xlabel('\omega_x');
    ylabel('\omega_y');
    title('FFT: Power Spectrum (log_{10})');
    colorbar;
    colormap jet;
    shading interp;
    axis tight;
    
    R = RAPS(Kx,Ky,log10(amplitude.^2));
    subplot(2,2,3);
    scatter(R(:,1),R(:,2));
    xlabel('\omega_x');
    ylabel('\omega_y');
    title('FFT: Radial Spectrum (log_{10})');
    axis tight;
end

function [result] = Fk(values)
    N = length(values);
    Ts = (values(N) - values(1)) / (N - 1);
    dF = 1 / (N * Ts);
    result = ((0:N-1) - ceil(N/2))*dF;
end

function [result] = RAPS(fx, fy, values)
    nx = length(fx);
    ny = length(fy);
    result = zeros(nx*ny,2);
    display(length(result));
    for j=1:ny,
        for i=1:nx,
            result((j-1)*nx + i,1) = (fy(i,j)^2 + fx(i,j)^2)^(0.5);
            result((j-1)*nx + i,2) = values(i,j);
       end
    end
end