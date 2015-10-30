function ft1d()
%
% function ft1d(datf)
%
% Perform 1D FFT on a 1D data profile and display
% resulting quantities for understanding the
% Fourier transform of different data types.
%
% Input:  datf: name of the data file (string variable)
%
% Requires the input data file in the following format:
%
%------------------------------------------------------
%
%N
%x(1),       y(1)
%x(2),       y(1)
%..
%x(N),       z(N)
%
%--------------------------------------------------------
%
%
%--------------------------------------------------------

% if nargin <1,
%    datf=input('  Name of data file: ');
% end;
datf=uigetfile('*.dat','Input 1D data file');

%====================================
% Step-0: Input the data
%====================================
fid = fopen(datf,'rt');
nx = fscanf(fid,'%d',1);
for i=1:nx,
   x(i) = fscanf(fid,'%f',1);
   y(i) = fscanf(fid,'%E',1);
end
fclose(fid);

%====================================
% Step-2: Expand the grid to the nearest
%         power of 2 
%====================================
%n2x = pow2(ceil(log2(nx)));
%n2y = pow2(ceil(log2(ny)));
n2x=nx;
%
%
% calculate angular frequencies: 
% Note: the frequencies are ordered from 0 to kx (Nyquist),
%       then from -kx+1 to -1
%       dxo and dyo are the frequency intervals
%
% Nyquist indices
kx = n2x/2;

% frequencies
dx=x(2)-x(1);
dxo=2*pi/(dx*n2x);
for ii=-kx+1:kx,
   omegax(ii+kx)=ii*dxo;
end;
%
%===============================
% Step-2: apply FFT 
%===============================
tmp=foldfft(y);
zft=fft(tmp);
zft=unfoldfft(zft);
%
amplitude=abs(zft);
real_part=real(zft);
imag_part=imag(zft);
phase=angle(zft);

%================================
% plot the data:
%=================================
figure('position', [0, 0, 600, 450]);
subplot(2,2,1),
plot(x, y);
axis tight;
xlabel('X (m)')
ylabel('f(x)')
title('Space-domain Data')

subplot(2,2,2),
semilogy(omegax,amplitude.^2);
axis tight;
xlabel('\omega (rad/m)')
ylabel('P(\omega)')
title('Power Spectrum')

subplot(2,2,3),
plot(omegax, real_part);
axis([omegax(1) omegax(nx) min(min(real_part),min(imag_part)) max(amplitude) ]);
xlabel('\omega (rad/m)')
ylabel('Real part')
title('Real part of FT')


subplot(2,2,4),
plot(omegax, imag_part);
axis([omegax(1) omegax(nx) min(min(real_part),min(imag_part)) max(amplitude) ]);
xlabel('\omega (rad/m)')
ylabel('Imaginary part')
title('Imaginary part of FT')

