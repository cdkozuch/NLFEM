function [b2,tau2,qBar2,alpha2,C2] = returnMap(f2,F2,b1,qBar1,alpha1,matProps)
%%RETURNMAP performs return mapping routine

kappa = matProps(1); C1 = matProps(3);
K0 = matProps(4); K1 = matProps(5); H = matProps(6);

fbar = (det(f2)^(-1/3))*f2;
bTrial = fbar*b1*fbar';

J = det(F2);
Up  = 0.5*kappa*(J - 1/J); %derivative of volumetric energy
CTrial = matStiffTen3D(F2,J,kappa,C1);

sTrial = 2*C1*deviatoric(bTrial);
xiTrial = sTrial - qBar1; %trial relative stress 
xiTrialMag = norm(xiTrial,'fro');

fTrial = xiTrialMag - sqrt(2/3)*(K0 + K1*alpha1); %yield function
if fTrial <= 0
    %still in elastic zone    
    tau2 = J*Up*eye(3) + sTrial;
    b2 = bTrial;
    qBar2 = qBar1;
    alpha2 = alpha1;  
    C2 = CTrial;
    return
end

muBar = (1/3)*(2*C1 + H)*trace(bTrial);
dGamma = (fTrial/2)/(muBar + K1/3);

nMat = xiTrial/xiTrialMag;
Lvb = (2/3)*dGamma*trace(bTrial)*nMat;

alpha2 = alpha1 + sqrt(2/3)*dGamma;
qBar2 = qBar1 + H*Lvb;
s = sTrial - 2*C1*Lvb;
tau2 = J*Up*eye(3) + s;
b2 = bTrial - Lvb;

beta0 = 1 + K1/(3*muBar);
beta1 = 2*muBar*dGamma/xiTrialMag;
beta2 = (2/3)*(1 - 1/beta0)*xiTrialMag*dGamma/muBar;
beta3 = 1/beta0 - beta1 + beta2;
beta4 = (1/beta0 - beta1)*xiTrialMag/muBar;

nMat = xiTrial/xiTrialMag;
nVec = reduceOrder21(nMat);
IVec = [ones(3,1);zeros(3,1)]; %eye(3) in vector form
tmp = nVec*reduceOrder21(deviatoric(nMat*nMat))';
symTmp = 0.5*(tmp + tmp');
matIden4 = diag([ones(1,3),0.5*ones(1,3)]); %4-tensor iden. as 2-tensor
II = [ones(3),zeros(3);zeros(3),zeros(3)]; %dyadic product of 2 eye(3)
CBar = 2*muBar*(matIden4 - (1/3)*II) - (2/3)*xiTrialMag*(nVec*IVec' + IVec*nVec');
C2 = CTrial - beta1*CBar - 2*muBar*beta3*(nVec*nVec') - 2*muBar*beta4*symTmp;

end


