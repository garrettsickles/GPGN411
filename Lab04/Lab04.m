function Lab04()
    filename = uigetfile('*.dat','Input 2D data file');
    file = fopen(filename, 'rt');
    ny = fscanf(file,'%d',1);
    nx = fscanf(file,'%d',1);
    x = zeros(nx);
    y = zeros(nx);
    z = zeros(nx);
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
    
    figure('position', [0, 0, 600, 750]);
    subplot(3,2,1);
    pcolor(y, x, z);
    xlabel('x (m)');
    ylabel('y (m)');
    title('Data'); 
    ylabel(colorbar, 'Amplitude')
    colormap jet;
    shading interp;
    axis tight;

    subplot(3,2,2);
    pcolor(Ky, Kx, log(amplitude.^2));
    xlabel('\omega_x');
    ylabel('\omega_y');
    title('FFT: Power Spectrum');
    ylabel(colorbar, 'log_{10}(Amplitude^{2})')
    colormap jet;
    shading interp;
    axis tight;
    
    % Radially averaged power spectrum
    R = RAPS(Kx,Ky,amplitude);
    
    % Linear interpolation
    x = min(R(:,1)):0.0001:max(R(:,1));
    y = interp1(R(:,1), R(:,2), x);
    
    % Trapezoidal integration
    ing = TrapInt(x, y, 200);
    
    subplot(3,2,3:6);
    scatter(R(:,1),log(R(:,2).^2),'.');
    hold on;
    plot(x, log(y.^2), 'color', [0.0,0.0,0.0]+0.7);
    hold on;
    plot(ing(:,1), log(ing(:,2).^2), 'color', 'red','marker','o','markerfacecolor','r');
    
    xlabel('\omega_r');
    ylabel('ln(Amplitude^2)');
    title(['Radially Averaged Power Spectrum for ', filename]);
    legend('Data','Interpolated','Integrated');
    axis tight;
end

% Produce a frequency spectrum
function [result] = Fk(values)
    N = length(values);
    Ts = (values(N) - values(1)) / (N - 1);
    dF = 1 / (N * Ts);
    result = ((0:N-1) - ceil(N/2))*dF;
end

% Radially Averaged Power Spectrum
function [avg] = RAPS(fx, fy, values)

    % Unwrap the 2D meshgrid into a list of pairs
    nx = length(fx);
    ny = length(fy);
    result = zeros(nx*ny,2);
    avg = zeros(2, 1);
    for j = 1:ny
        for i = 1:nx
            result((j-1)*nx + i,1) = (fy(i,j)^2 + fx(i,j)^2)^(0.5);
            result((j-1)*nx + i,2) = values(i,j);
       end
    end
    
    % Sort the result 
    result = sortrows(result,1);
    nr = length(result);
    
    % Base Cases
    first = 1;
    last = 1;
    sum = result(first,2);
    unique = 1;
    
    % N order data shrinking and averaging
    while last < nr
        last = last + 1;
        % If next point has same radial frequency
        if result(first, 1) == result(last, 1)
            % Add it to the running sum
            sum = sum + result(last, 2);
        else
            % Set the next unique frequency and amplitude
            avg(unique, 1) = result(first, 1);
            avg(unique, 2) = sum / (last - first);
            
            % Set first and last to the same
            % Reset sum to the next value
            first = last;
            sum = result(last, 2);
            unique = unique + 1;
        end
    end
end

% Trapezoidal integration
function result = TrapInt(x, y, spacing)
    py = 0;
    px = 0;
    index = 1;
    first = 1;
    result = zeros(floor(length(x)/spacing), 2);
    for i = 1:(length(x)-1)
        py = py + (0.5)*(x(i+1)-x(i))*(y(i)+y(i+1));
        px = px + x(i);
        if mod(i, spacing) == 0
            result(index, 1) = px / (spacing);
            result(index, 2) = py / (x(i+1) - x(first));
            first = i + 1;
            px = 0;
            py = 0;
            index = index + 1;
        end
    end
end