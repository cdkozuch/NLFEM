function B = reduceOrder21(A)
%%REDUCEORDER21 converts symmetric 2-tensor to vector notation
B = [A(1,1) A(2,2) A(3,3) A(1,2) A(2,3) A(1,3)]';
end