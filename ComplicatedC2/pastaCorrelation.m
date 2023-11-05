function [q,C2T] = pastaCorrelation
a=110;
b=-26;
c=24;
big_lambda=1.25;
small_lambda=12.35;  % not a 100% exact
c_p = 0.8;
eV=1.602176634*10^(-19);
eps0=8.8541878128*10^(-12);
alpha_const=(1.602176634*10^(-19))^2/(4*pi*eps0)*1*10^12*(1*10^(-6)/eV); % in MeV*fm, EXTREMELY SMALL, SO SMALL THAT V_NN=V_PP ???
beta=(1)^(-1); % temperature in eV, in fact imply that it is beta. So, beta= 1/1 MeV^-1
n= 0.05; % in 1/fm^3   

m=10000;
delr = 0.001;
r = (-m:m)*delr;
rpos = r(r>=0);

V_pp=a*exp(-(rpos).^2/(big_lambda))+(b+c)*exp(-(rpos).^2/(2*big_lambda))+alpha_const./abs(rpos).*exp(-abs(rpos)/small_lambda);
V_np=a*exp(-(rpos).^2/(big_lambda))+(b-c)*exp(-(rpos).^2/(2*big_lambda));
V_nn=a*exp(-(rpos).^2/(big_lambda))+(b+c)*exp(-(rpos).^2/(2*big_lambda));

C2pp = (n*(exp(-beta*V_pp)-1));
C2np = (n*(exp(-beta*V_np)-1));
C2nn = (n*(exp(-beta*V_nn)-1));

[~,C2ppk] = C2_radial_transform(C2pp,r,rpos,delr,m);
[~,C2npk] = C2_radial_transform(C2np,r,rpos,delr,m);
[q,C2nnk] = C2_radial_transform(C2nn,r,rpos,delr,m);

C2T = (c_p^2*C2ppk+(1-c_p)^2*C2nnk+2*c_p*(1-c_p)*C2npk);

% Remove everything prior to the first trough
[~,qmin] = min(C2T);
C2T = C2T(q>q(qmin));
q = q(q>q(qmin));
q=q-q(1);
end