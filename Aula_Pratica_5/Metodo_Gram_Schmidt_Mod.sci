function [Q,R]=qr_GSM(A)

[m,n]=size(A);
Q = eye(m,n);  //ortogonal (Q^TQ=I)
R = eye(n,n);  //triangular superior

for j=1:n
    v = A(1:m,j);
    for i=1:j-1
        R(i,j) = Q(1:m,i)' * v; //subtrai a projeção da última versão de v
        v = v - R(i,j) * Q(1:m,i); //sobre q_1, ..., q_j-1
    end
    R(j,j) = norm(v);
    Q(1:m,j) = v / R(j,j);
    
end

endfunction
