function [KT,Fint,F2,b2,sigma2,qBar2,alpha2,C2] = globTanStiff(elements,coors,U1,U2,F1,b1,qBar1,alpha1,matProps)
%%GLOBTANSTIFF computes global tangent stiffness and residual force

nten = 6; %number of unique stress/strain values in 3D

ndofn = size(coors,2); %dofs per node
nnode = size(coors,1); %number of nodes
nnel = size(elements,2); %nodes per elements
nel = size(elements,1); %number of elements
ndoft = nnode*ndofn; %total number of degrees of freedom
ndofel = nnel*ndofn; %number of degrees of freedom per element

gp = (1/sqrt(3))*[-1 1]; %gauss points
ngp = numel(gp); %number of guass points
nintpel = ngp^ndofn; %number of integration points per element
nintpt = nintpel*nel; %number of total integration points

%preallocate return variables
F2 = zeros(nten,nintpt);
b2 = zeros(nten,nintpt);
sigma2 = zeros(nten,nintpt);
qBar2 = zeros(nten,nintpt);
alpha2 = zeros(1,nintpt);
C2 = zeros(nten,nten,nintpt);

KT = spalloc(ndoft,ndoft,ndoft); %global stiffness matrix
Fint = spalloc(ndoft,1,ndoft); %residual force vector

%create global stiffness and tangent stiffness matrices
intpi = 1; %indegration point index
for iter=1:nel
    ind = globElemInd(elements,iter,ndofn); %global indices for element
    u = reshape(U1(ind),3,8); %displacement
    uOld = reshape(U2(ind),3,8);
    du = u - uOld;
    
    refElemCoors = reshape(coors(:,:,iter),8,3);
    elemCoors = refElemCoors + u';
    
    ke = zeros(ndofel); %element stiffness matrix
    fe = zeros(ndofel,1); %element internal force vector
    
    for i=1:ngp
      xi = gp(i); %first natural coordinate
      
      for j=1:ngp
        eta = gp(j); %second natural coordinate
        
        for k=1:ngp
          zeta = gp(k); %third natural coordinate

          natcoors = [xi,eta,zeta];            
          [B,B2,DN,J] = bmatHex8(elemCoors,natcoors);
          detJ = det(J);
          f2 = eye(3) + du*DN';
          
          F = increaseOrder12(F1(:,intpi));
          F = f2*F;
          b = increaseOrder12(b1(:,intpi));
          qBar = increaseOrder12(qBar1(:,intpi));
          alpha = alpha1(intpi);
          
          [b,tau,qBar,alpha,C] = returnMap(f2,F,b,qBar,alpha,matProps);
          sigmaMat = tau/det(F);

          PK2 = inv(F)*tau*inv(F'); %2nd piola kirchhoff stress
          PK2p = kron(eye(3),PK2);
            
          ke = ke + detJ*(B'*C*B + B2'*PK2p*B2); %element stiffness matrix
          fe = fe + detJ*B'*reduceOrder21(PK2); %element force vector
          
          F2(:,intpi) = reduceOrder21(F);
          b2(:,intpi) = reduceOrder21(b);
          sigma2(:,intpi) = reduceOrder21(sigmaMat);
          qBar2(:,intpi) = reduceOrder21(qBar);
          alpha2(intpi) = alpha;
          C2(:,:,intpi) = C;
          
          intpi = intpi + 1;
        end
      end
    end

    Fint(ind) = Fint(ind) + fe; %update global residual force vector
    KT(ind,ind) = KT(ind,ind) + ke; %update global stiffness matrix
end

end