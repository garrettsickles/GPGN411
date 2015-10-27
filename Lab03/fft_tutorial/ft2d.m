function ft2d()
%
% function ft2d()
%
% Perform 2D FFT on a 2D data map and display
% resulting quantities for understanding the
% Fourier transform of different data types.
%
% Input:  datf: name of the data file (string variable)
%
% Requires the input data file in the following format:
%
%------------------------------------------------------
%
%N_east, N_north
%y(1), x(1),       z(1,1)
%y(1), x(2),       z(1,2)
%..
%y(1), x(N_north), z(1, N_north)
%y(2), x(1),       z(2,1)
%...
%y(N_east), x(N_north), z(N_east, N_north)
%--------------------------------------------------------
%
%
%--------------------------------------------------------


datf=uigetfile('*.dat','Input 2D data file');


%====================================
% Step-0: Input the data
%====================================
fid = fopen(datf,'rt');
ny = fscanf(fid,'%d',1);
nx = fscanf(fid,'%d',1);
for j=1:ny,
	for i=1:nx,
		y(i,j) = fscanf(fid,'%f',1);
		x(i,j) = fscanf(fid,'%f',1);
        z(i,j) = fscanf(fid,'%E',1);
   end
end
fclose(fid);

%====================================
% Step-2: Expand the grid to the nearest
%         power of 2 
%====================================
%n2x = pow2(ceil(log2(nx)));
%n2y = pow2(ceil(log2(ny)));
n2x=nx;
n2y=ny;
%
%
% Padding with zero:
% Create a zero matrix of the appropriate dimensions
% and fill the upper-left portion with the input data
%
%lx1=fix((n2x-nx)/2);
%lx2=n2x-nx-lx1;
%ly1=fix((n2y-ny)/2);
%ly2=n2y-ny-ly1;

%zexpand=zeros(n2x,n2y);
%zexpand(lx1+1:lx1+nx,ly1+1:ly1+ny)=z;
%
% calculate angular frequencies: 
% Note: the frequencies are ordered from 0 to kx (Nyquist),
%       then from -kx+1 to -1
%       dxo and dyo are the frequency intervals
%
% Nyquist indices
kx = n2x/2;
ky = n2y/2;

% frequencies
dx=x(2,2)-x(1,1);
dy=y(2,2)-y(1,1);
dxo=2*pi/(dx*n2x);
dyo=2*pi/(dy*n2y);
for ii=-kx+1:kx,
   omegax(ii+kx)=ii*dxo;
end;
for ii=-ky+1:ky,
   omegay(ii+ky)=ii*dyo;
end;
[Omegay,Omegax]=meshgrid(omegay,omegax);
%
%===============================
% Step-2: apply FFT 
%===============================
tmp=foldfft(z);
zft=fft2(tmp);
zft=unfoldfft(zft);
%
amplitude=abs(zft);
real_part=real(zft);
imag_part=imag(zft);
phase=angle(zft);

%================================
% plot the data:
%=================================
xmin=min(min(x));
xmax=max(max(x));
ymin=min(min(y));
ymax=max(max(y));

pmax=max(omegax);
pmin=min(omegax);
qmax=max(omegay);
qmin=min(omegay);

tmax=max(max(amplitude));
tminr=min(min(real_part));
tmini=min(min(real_part));
tmin=min(tminr,tmini);

figure('position', [0, 0, 600, 450]);
subplot(2,2,1),
pcolor(y, x, z);
xlabel('x (m)');
ylabel('y (m)');
title('Data');
colorbar;
colormap jet;
shading interp;
axis tight;

subplot(2,2,2),
pcolor(Omegay, Omegax, log10(amplitude.^2));
xlabel('\omega_x');
ylabel('\omega_y');
title('FFT: Power Spectrum (log10)');
amax=2*log10(tmax);
amin=amax-10.0;
caxis([amin amax]);
colorbar;
colormap jet;
shading interp;
axis tight;

subplot(2,2,3),
pcolor(Omegay, Omegax, real_part);
xlabel('\omega_x');
ylabel('\omega_y');
title('FFT: Real Part');
shading interp;
axis tight;
colormap jet;
caxis([tmin, tmax]);
colorbar;

subplot(2,2,4),
pcolor(Omegay, Omegax, imag_part);
xlabel('\omega_x');
ylabel('\omega_y');
title('FFT: Imaginary Part');
shading interp;
axis tight;
colormap jet;colormap jet;
caxis([tmin, tmax]);
colorbar;

