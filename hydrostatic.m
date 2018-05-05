function B = hydrostatic(A)
%%HYDROSTATIC returns hydrostatic part of A
B = (1/3)*trace(A)*eye(3);
end