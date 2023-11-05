

function [ operators_h, operators_d ] = prepareOperators( params, gpu, dx, dy, dz)

%params=createParams();
%PREPAREOPERATORS
%params=createParams(); %even if we want to use previously made params, we need to use the params generator script for the location of data output


L=params.L;
%dx=params.dx;
dt=params.dt;
%Open the txt file containing the k-array and the effective correlation
%arrays and transfer it into a array
% File1=fopen(params.name_k_file);
% File2=fopen(params.name_effectivec2_file);
% k_correlation= fscanf(File1, "%f");
% C_effective=fscanf(File2,"%f");
% fclose(File1);
% fclose(File2);

%prepare the k array for X coordinate
fac_x = 2*pi / (dx*L);
fac_y = 2*pi / (dy *L);
fac_z = 2*pi / (dz*L);

kx_list = ifftshift(fac_x*gpuArray(linspace(-L/2,L/2-1,L)));
ky_list = ifftshift(fac_y*gpuArray(linspace(-L/2,L/2-1,L)));
kz_list = ifftshift(fac_z*gpuArray(linspace(-L/2,L/2-1,L)));

%disp(k_list)

if params.is3D
    [ki,kj,kl] = meshgrid(kx_list,ky_list,kz_list);
else
    [ki,kj] = meshgrid(kx_list,ky_list);
end



%prepare the linear and non-linear operators

if params.is3D
    kLap = -(ki.*ki + kj.*kj + kl.*kl);
else
    kLap = -(ki.*ki + kj.*kj);
end


%Return the effective correlation function value for the point obtain by

% If we want to add a strain flag, put it here.


[k_correlation,C_effective] = pastaCorrelation();
C_2_spline=spline(k_correlation,C_effective);
C_2 = gpuArray(ppval(C_2_spline,sqrt(-kLap)));


% if max(k_correlation,[], "all")<max(-kLap,[], "all")
%     disp('Max of simulation:')
%     disp(max(-kLap,[], "all"))
%     disp('Max of file')
%     disp(max(k_correlation,[], "all"))
%     error('DANGER, the k range given in the file is too small compare to the k range that you are trying to simulate')
% end


opCk =1-C_2;

opL = exp(kLap.*opCk*dt);

ind = (opCk==0);

opN =   ((exp(kLap.*opCk*dt) - 1.) ./ opCk);
opN(ind) = kLap(ind)*dt;

if params.dispOutput
    plot(1:L/2,opCk(1:L/2,1,1),1:L/2,opL(1:L/2,1,1),1:L/2,opN(1:L/2,1,1));
end


operators_h = struct();
operators_d = struct();

operators_h.kLap = kLap;
operators_h.opL = opL;
operators_h.opN = opN;

operators_d.opL = gpuArray(opL);wait(gpu);
operators_d.opN = gpuArray(opN);wait(gpu);


end

