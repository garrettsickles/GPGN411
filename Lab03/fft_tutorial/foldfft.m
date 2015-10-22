function aout=foldfft(ain)
%
% function aout=fold(ain)
%
% fold Fourier transform to shift
% the DC component to lower-left corner
%
nd = ndims(ain);
idx = cell(1, nd);
for k = 1:nd
    nx = size(ain, k);
    kx = floor(nx/2);
    if kx>1
       idx{k} = [kx:nx 1:kx-1];
    else
       idx{k} = [1];
    end
end
aout = ain(idx{:});