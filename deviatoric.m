function B = deviatoric(A)
%%DEVIATORIC returns deviatoric part of A
B = A - hydrostatic(A);
end

