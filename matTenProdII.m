function B = matTenProdII()
%%MATIDEN4 2-tensor representation of tensor product of eye(3) and eye(3)

A = zeros(3,3,3,3);
for i=1:3
    for j=1:3
        dij = i==j;
        for k=1:3
            for l=1:3
                A(i,j,k,l) = dij*(k==l);
            end
        end
    end
end
B = reduceOrder42(A);

end