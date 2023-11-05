function [ fields_d, status ] = stepPFC( params, gpu, fields_d, operators_d, status )
%STEPPFC 

psiN = -1/2*(fields_d.psi.^2) + 1/3*(fields_d.psi.^3); wait(gpu);

psiN = fftn(psiN);wait(gpu);


fields_d.psi_F = operators_d.opL.*fields_d.psi_F + operators_d.opN.*psiN;wait(gpu);

fields_d.psi=real(ifftn(fields_d.psi_F));wait(gpu);

end

