%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script executes the non-linear finite element method for a mesh 
% composed of one 8-node hexahedral element with combined isotropic and
% kinematic hardening. The material model is neo-hookean. The element is 
% placed under shear by applying a fixed dispalacement corresponding to 2% 
% strain in the e1 direction to nodes 3, 4, 7, and 8 while pinning nodes 1, 
% 2, 5, and 6.
%
%      5-------8
%    / |     / |
%   6--|----7  |
%   |  1----|--4
%   | /     | /
%   2-------3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% Get nodal coordinates from data files
load 'nodes.dat';
load 'elements.dat';
nnode = size(nodes,1); %number of nodes 
nel = size(elements,1); %number of elements
nnel = size(elements,2); %number of nodes per element
ndofn = size(nodes,2); %number of degrees of freedom per node
coors = zeros(nnel,ndofn,nel); %nodal coordinates for each element
ndoft = nnode*ndofn; %total number of degrees of freedom
ndofel = nnel*ndofn; %number of degrees of freedom per element

for i=1:nel
  for j=1:nnel
    coors(j,:,i) = nodes(elements(i,j),:);
  end
end

%% elastic material properties
E = 1.199; %elastic modulus (MPa)
mu = 0.412; %shear modulus (MPa)
lambda = mu*(E - 2*mu)/(3*mu - E); %first lame constant
kappa = lambda + (2/3)*mu; %bulk modulus

%% hyperelastic material properties
C1 = 0.165; %neo-hookean constant (MPa)

%% Plastic material properties
K0 = 0.03; %yield stress
K1 = 0.5; %isotropic hardening modulus
H = 0.01; %kinematic hardening modulus
matProps = [kappa,mu,C1,K0,K1,H];
      
%% Define boundary conditions
utotal = 0.2; %total displacement of top surface in loading direction
nstep1 = 100; %number of steps until displacement total
nstep2 = 100; %final number of steps
du = utotal/nstep1; %displacement increment

fixedDofs = [1 2 3 4 5 6 13 14 15 16 17 18 ...
             7 10 19 22]; %fixed displacement DoFs
dFixedVals = [zeros(1,12) du*ones(1,4)]; %fixed displacment increment values

%fixedDofs = [1 2 3 10 12 13 22 4 7 16 19]; %fixed displacement DoFs
%dFixedVals = [zeros(1,7) du*ones(1,4)]; %fixed displacment increment values

nFixedDofs = numel(fixedDofs);
freeDofs = setdiff(1:ndoft,fixedDofs);
Fext = spalloc(ndoft,1,ndoft); %no external forces

%% Initialize plastic history variables
nintpel = 8; %number of integration points per element (2x2x2 quadrature)
nintpt = nintpel*nel; %number of integration points
nten = 6; %number of unique strain/stress tensor elements
ep = zeros(nten,nintpt); %plastic strain
qbar = zeros(nten,nintpt); %kinematic hardening
alpha = zeros(nintpt); %effective plastic strain

%% Preallocate history storage
sigmaRec = cell(1); %stress history
sigmaRec{1} = zeros(nten,nintpt);
FRec = cell(1);
FRec{1} = repmat([ones(3,1);zeros(3,1)],1,nintpt);
bRec = cell(1);
bRec{1} = repmat([ones(3,1);zeros(3,1)],1,nintpt);
alphaRec = zeros(nstep2+1,nintpt); %effective plastic strain history
qBarRec = cell(1); %kinematic hardening history
qBarRec{1} = qbar;
Urec = zeros(ndoft,nstep2+1); %displacement (columns are time steps)    
Crec = cell(1); %material stiffness tensor history
C = matStiffTen3D(eye(3),1,kappa,C1);
Crec{1} = repmat(C,1,1,nintpt);

%% Step through analysis
maxni = 200; %maximum number of newton iterations
tolFac = 1e-12; %tolerance factor
killFlag = 0;

tic
for t=1:nstep2
    U2 = zeros(ndoft,1);
    U2(fixedDofs) = t*dFixedVals;
    
    U1 = Urec(:,t);
    F1 = FRec{t}; b1 = bRec{t};
    qBar1 = qBarRec{t}; alpha1 = alphaRec(t,:);
    
    %Newton iterations
    for ni=1:maxni
      %Build tangent stiffness and internal force
      [KT,Fint,F2,b2,sigma2,qBar2,alpha2,C2] = ...
          globTanStiff(elements,coors,U2,U1,F1,b1,qBar1,alpha1,matProps);     
      
      %Compute residual force
      RF = Fext - Fint;
      
      %Apply boundary conditions
      KT(fixedDofs,:) = zeros(nFixedDofs,ndoft);
      KT(fixedDofs,fixedDofs) = eye(nFixedDofs);

      %Evaluate condition
      res = norm(RF(freeDofs));
         
      if ni==1
          tol = tolFac*res;
      elseif res <= tol
        fprintf('Time step: %i;\tIteration: %i;\tResidual: %g\n',t,ni,res)
        break
      elseif ni==maxni || res > 1
        fprintf('Convergence failure!\n')
        nstep2 = t;
        killFlag = 1;
        break
      end

      if killFlag
          break
      end
      
      %Update displacement guess
      U2(freeDofs) = U2(freeDofs) + ...
        KT(freeDofs,freeDofs)\RF(freeDofs);
    end 

    %record history variables
    Urec(:,t+1) = U2;
    bRec{t+1} = b2;
    FRec{t+1} = F2;
    sigmaRec{t+1} = sigma2;
    Crec{t+1} = C2;
    qBarRec{t+1} = qBar2;
    alphaRec(t+1,:) = alpha2;
end

%% Plot results
intpi = 3; %integration point at which to plot

%plot stress vs displacement
d1 = Urec(10,1:nstep2+1); %displacement in e1 direction at node 4
s11 = zeros(1,nstep2+1); s22 = zeros(1,nstep2+1);
s33 = zeros(1,nstep2+1); s12 = zeros(1,nstep2+1);
for t=1:nstep2+1
  s11(t) = sigmaRec{t}(1,intpi); s22(t) = sigmaRec{t}(2,intpi);
  s33(t) = sigmaRec{t}(3,intpi); s12(t) = sigmaRec{t}(4,intpi);
end
figure
plot(d1,s11,d1,s22,d1,s33,d1,s12,'LineWidth',2)
legend('s11','s22','s33','s12','Location','NorthWest')
xlabel('displacement')
ylabel('sigma')

%plot stress vs time
tf1 = 1; dt = tf1/nstep1; t0 = 0; tf2 = dt*nstep2;
tvec = t0:dt:tf2; %time vector
figure
plot(tvec,s11,tvec,s22,tvec,s33,tvec,s12,'LineWidth',2)
legend('s11','s22','s33','s12','Location','NorthWest')
xlabel('time')
ylabel('sigma')

%plot alpha
figure
plot(tvec,alphaRec(1:nstep2+1,intpi),'LineWidth',2)
xlabel('time')
ylabel('alpha')

%plot qbar
q11 = zeros(1,nstep2+1); q22 = zeros(1,nstep2+1); 
q33 = zeros(1,nstep2+1); q12 = zeros(1,nstep2+1);
for t=1:nstep2+1
  q11(t) = qBarRec{t}(1,intpi); q22(t) = qBarRec{t}(2,intpi);
  q33(t) = qBarRec{t}(3,intpi); q12(t) = qBarRec{t}(4,intpi);
end
figure
plot(tvec,q11,tvec,q22,tvec,q33,tvec,q12,'LineWidth',2)
legend('qbar11','qbar22','qbar33','qbar12','Location','NorthWest')
xlabel('time')
ylabel('qbar')

%plot C
C1111 = zeros(1,nstep2+1); C2222 = zeros(1,nstep2+1); C1212 = zeros(1,nstep2+1);
for t=1:nstep2+1
  C1111(t) = Crec{t}(1,1,intpi);
  C2222(t) = Crec{t}(2,2,intpi);
  C1212(t) = Crec{t}(4,4,intpi);
end
figure
plot(tvec,C1111,tvec,C2222,tvec,C1212,'LineWidth',2)
legend('C1111','C2222','C1212')
xlabel('time')
ylabel('C')


