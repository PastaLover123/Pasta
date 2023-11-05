function [ operators_h, operators_d ] = prepareOperators( params, gpu, dx, dy, dz )
%PREPAREOPERATORS


L=params.L;

%dx=params.dx;
dt=params.dt;
Bx=params.Bx;
Bl=params.Bl;

%prepare the k array for X coordinate
facx = 2*pi / (dx*L);
facy = 2*pi / (dy*L);
facz = 2*pi / (dz*L);
for i=0:(L-1)
    if (i<L/2)
        kx_list(i+1) = i*facx;
        ky_list(i+1) = i*facy;
        kz_list(i+1) = i*facz;
    else
        kx_list(i+1) = (i - L)*facx;
        ky_list(i+1) = (i - L)*facy;
        kz_list(i+1) = (i - L)*facz;
    end
end

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

opCk = Bl + Bx*(2 * kLap + kLap.*kLap); %ken-style C2

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

