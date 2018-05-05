function D = matStiffTen3D(b,J,kappa,C1)
%%MATSTIFFTEN3D computes material stiffness tensor for neo-hookean material

mu = 2*C1;
muBar = (1/3)*mu*trace(b);
s = mu*deviatoric(b);
sMag = norm(s,'fro');
nMat = s/sMag;

matIden4 = diag([ones(1,3),0.5*ones(1,3)]); %4-tensor iden. as 2-tensor
II = [ones(3),zeros(3);zeros(3),zeros(3)]; %dyadic product of 2 eye(3)
IVec = reduceOrder21(eye(3)); %eye(3) in vector form
nVec = reduceOrder21(nMat);
Cbar = 2*muBar*(matIden4 - (1/3)*II) - (2/3)*sMag*(nVec*IVec' + IVec*nVec');
D = kappa*(J^2)*II - kappa*(J^2 - 1)*matIden4 + Cbar;
end