function [x] = Resolve_Lx(L,b)

[t]=size(L,1);
x=zeros(t,1);

// Calcula x, sendo Lx=b

x(1)=b(1)/L(1,1);
for i=2:t
    x(i)=(b(i)-L(i,1:i)*x(1:i))/L(i,i);
end

endfunction

function [xk,N,k,rk]=Gauss_Seidel_Method_2(A,b,x0,E,M,norma)
   
k = 0
xk=x0     
rk=norm(b-A*xk,norma) //Resíduo da operação    
D=diag(diag(A)) //Diagonal de A
LPU = A-D //L plus U
L=tril(A,-1)
U=triu(A,1)


while k<M | N>E
    xk1=Resolve_Lx(L+D,b-U*xk)
    N=norm(xk1-xk,norma) //norma da diferença entre as aproximações
    xk=xk1
    k=k+1
    if k>M
        break
    end    
end

endfunction
