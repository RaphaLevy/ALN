function [Q,R,P]=qr_GSP(A)

[m,n] = size(A);
Q = eye(m,n);  //ortogonal (Q^TQ=I)
R = eye(n,n);  //triangular superior
P = eye(n,n);  //matriz de pivoteamento
X = A

aux = [];                             //Pega a coluna de maior norma
for i=1:n
    aux(i) = norm(A(:,i));
end

[x, pos]=max(norm(A(:,i)));
v = A(:,pos);
R(1,1) = norm(v);
Q(:,1) = v / norm(v);

if pos~=1 then                         //Troca de colunas
    A(:,[1, pos])=A(:,[pos, 1]);
    P(:,[1, pos])=P(:,[pos, 1]);
end

for j=2:n
    A(:,j) = A(:,j) - ((Q(:,j-1)' * A(:,j)) / (Q(:,j-1)' * Q(:,j-1))) * Q(:,j-1);
        for l=2:n
        maximo = max(norm(A(:,l)));
    end    
    a=(j-1) + maximo;
    A(:,[j, a])=A(:,[a, j]); //Erro de índice inválido?
    P(:,[j, a])=P(:,[a, j]);
    
    v = A(:,j);
    for i=1:j-1
        qTa_ij = Q(:,i)' * v;
        k = [];
        for z=1:n
            k(z) = Q(:,i)' * X(:,z);
        end
        R(i,j) = k' * P(:,j);
        v = v - qTa_ij * Q(:,i);
    end
    
    Q(:,j)=v/norm(v,2);
    R(j,j) = norm(v);
end

endfunction
