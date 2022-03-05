function [xk,N,k,rk]=Jacobi_Method(A,b,x0,E,M,norma)
       
k = 0
xk=x0         
rk=norm(b-A*xk,norma) //Resíduo da operação    
D=diag(diag(A)) //Diagonal de A
LPU = A-D //L plus U
invD=diag(1./diag(D))  //Inversa da matriz D
MJ=(-1)*invD*LPU //Matriz de Jacobi
cJ=invD*b //Vetor de Jacobi 

while k<M | N>E
    xk1=MJ*xk+cJ
    N=norm(xk1-xk,norma) //norma da diferença entre as aproximações
    xk=xk1
    k=k+1
    if k>M
        break
    end  
end

endfunction
