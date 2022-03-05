function [U,R]=qr_House(A)

[m,n]=size(A);
U = zeros(m,n);
R = eye(m,n);   //triangular superior

for k=1:n
    x = A(k:m,k);
    if length(x) > 1 then
        if x(1) < 0 then              //calcula x-Hx
           x(1) = x(1) - norm(x);     //calcula x-Hx
        else                          //calcula x-Hx         
           x(1) = x(1) + norm(x);     //calcula x-Hx  
        end
        u = x / norm(x);
        U(k:m,k) = u;
        A(k:m,k:n) = A(k:m,k:n) - 2 * u * (u' * A(k:m,k:n)); //Passo da triangularização
    else
        U(k:m,k) = 1;
        A(k:m,k:n) = x;
    end
end

R = triu(A);

endfunction

function [Q]=constroi_Q_House(U)
 
[m,n]=size(U);   
I=eye(m,m);
u_1 = U(:,1);
Q = I - 2 * u_1 * u_1';

for x=2:n
    u_x = U(:,x);
    Q_x = I - 2 * u_x * u_x';
    Q = Q * Q_x;
end

endfunction


