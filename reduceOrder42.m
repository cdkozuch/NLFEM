function B = reduceOrder42(A)
%%REDUCEORDER21 converts symmetric 4-tensor to 2-tensor notation
B = [ A(1,1,1,1) A(1,1,2,2) A(1,1,3,3) A(1,1,1,2) A(1,1,2,3) A(1,1,1,3);
      A(2,2,1,1) A(2,2,2,2) A(2,2,3,3) A(2,2,1,2) A(2,2,2,3) A(2,2,1,3);
      A(3,3,1,1) A(3,3,2,2) A(3,3,3,3) A(3,3,1,2) A(3,3,2,3) A(3,3,1,3);
      A(1,2,1,1) A(1,2,2,2) A(1,2,3,3) A(1,2,1,2) A(1,2,2,3) A(1,2,1,3);
      A(2,3,1,1) A(2,3,2,2) A(2,3,3,3) A(2,3,1,2) A(2,3,2,3) A(2,3,1,3);
      A(1,3,1,1) A(1,3,2,2) A(1,3,3,3) A(1,3,1,2) A(1,3,2,3) A(1,3,1,3) ];
end