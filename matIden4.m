function B = matIden4()
%%MATIDEN4 2-tensor representation of 4-tensor identity

A = zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            dik = i==k;
            djk = j==k;
            for l=1:3
                A(i,j,k,l) = 0.5*(dik*(j==l) + (i==l)*djk);
            end
        end
    end
end
B = reduceOrder42(A);

end
