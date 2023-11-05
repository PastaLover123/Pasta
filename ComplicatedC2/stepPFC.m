
function [ fields_d, status ] = stepPFC( params, gpu, fields_d, operators_d, status )
%STEPPFC 
noise =1;
if noise == 1
    X = params.sigma*randn(params.L,params.L,params.L);
    Y = params.sigma*randn(params.L,params.L,params.L);
    Z = params.sigma*randn(params.L,params.L,params.L);
    R = gpuArray((X - circshift(X,[-1,0,0]) + Y - circshift(Y,[0,-1,0]) + Z - circshift(Z,[0,0,-1]))/params.dx);
    Rk = fft2(R);
end

psiN = -1/2*(fields_d.psi.^2) + 1/3*(fields_d.psi.^3)+R; wait(gpu); % extra term not dep on B

psiN = fftn(psiN);wait(gpu);

%disp("Before:")
%disp('Max F:')
%disp(max(abs(fields_d.psi_F),[],"all"))
fields_d.psi_F = operators_d.opL.*fields_d.psi_F + operators_d.opN.*psiN;wait(gpu);

fields_d.psi=real(ifftn(fields_d.psi_F));wait(gpu);
%disp('After')
%disp('Max real:')
%disp(max(abs(fields_d.psi),[],"all"))
%disp('Max F:')
%disp(max(abs(fields_d.psi_F),[],"all"))
%disp('max operator L:')
%disp(max(abs(operators_d.opL),[],"all"))
%disp('max operator N:')
%disp(max(abs(operators_d.opN),[],"all"))


%N is the non linear operator calculated at each time step
%L is the linear operator pre-calculated at the start of the simulation


end

