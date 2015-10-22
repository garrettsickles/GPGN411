function aout=unfoldfft(ain)
%
% function aout=fold(ain)
%
% unfold Fourier transform to shift
% the DC component to centre
%
nd = ndims(ain);
idx = cell(1, nd);
for k = 1:nd
    nx = size(ain, k);
    kx = floor(nx/2);
    if kx>1
       idx{k}=[kx+2:nx 1:kx+1];
    else
       idx{k}=[1];
    end
end
aout = ain(idx{:});