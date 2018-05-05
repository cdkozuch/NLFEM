function B = increaseOrder12(A)
%%INCREASEORDER12 converts vector notation to symmetric 2-tensor
B = [ A(1) A(4) A(6);
      A(4) A(2) A(5);
      A(6) A(5) A(3) ];
end
     