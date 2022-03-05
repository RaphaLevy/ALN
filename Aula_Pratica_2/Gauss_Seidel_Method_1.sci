function [xk,N,k,rk]=Gauss_Seidel_Method_1(A,b,x0,E,M,norma)
         
k = 0
xk=x0         
rk=norm(b-A*xk,norma) //Resíduo da operação    
D=diag(diag(A)) //Diagonal de A
LPU = A-D //L plus U
invD=diag(1./diag(D))  //Inversa da matriz D
L=tril(A,-1)
U=triu(A,1)
MG=(-1)*(inv(L+D))*(U) //Matriz de Jacobi
cG=(inv(L+D))*(b) //Vetor de Jacobi

while k<M | N>E
    xk1=MG*xk+cG
    N=norm(xk1-xk,norma) //norma da diferença entre as aproximações
    xk=xk1
    k=k+1
    if k>M
        break
    end    
end

endfunction

function [DiagDom]=Diagonal_Dominante(A)
    
[n]=size(A,1);    
    
for i=1:n
    for j=1:n
        if abs(A(i,i))<sum(abs(A(i,j)))
                     break
        else
            return DiagDom             
        end
    end
end    
endfunction
