function [q,C2q] = C2_radial_transform(Frpos,r,rpos,delr,m)
rFrpos = rpos.*Frpos;
% create antisymmetric function with r=0 at the center
% see plot 2
rFr = [-fliplr(rFrpos(2:end)) rFrpos];
N = length(r);
qFq = -imag(fftshift(fft(ifftshift(rFr))))*delr;
delq = 2*pi/(N*delr);
q = (-m:m)*delq;
C2q = qFq((m+2):end)./q((m+2):end);
q=q((m+2):end);
%C2q_Inverse = (1/(2*pi))*N*imag(fftshift(ifft(ifftshift(qFq))))*delq./r;
end