function [lambda,x1,k,n_erro]=Metodo_potencia_1(A,x0,epsilon,M)
    
//VARIÁVEIS DE ENTRADA    
//A: matriz nxn diagonalizável com autovalor dominate lambda
//x0: vetor NÃO NULO a ser usado como aproximação inicial do autovetor dominante    
//epsilon: precisão a ser usada no critério de parada
//M: número máximo de iterações

//VARIÁVEIS DE SAÍDA
//lambda: autovalor dominante de A
//x1: autovetor unitário (norma infinito) correspondente a lambda
//k: número de iterações necessárias para a convergência
//n_erro: norma infinito do erro
    
k=1
[m,i]=max(abs(x0))
x0=x0/x0(i) //coord maior modulo x0  
x1=A*x0 //aproximacao do autovetor dominante
while k<=M
    [m2,i2]=max(abs(x1)) //coord maior modulo x1; aproximacao do autovalor dominante
    lambda=x1(i2)
    x1=x1/lambda
    n_erro=norm(x1-x0,%inf)
    if n_erro < epsilon
        disp("O metodo convergiu com a precisao passada")
        break
    end   
    x0=x1
    x1=A*x0    
    k=k+1
end 

if k>M    
    disp("O numero de iteracoes passou do maximo de iteracoes permitido")   
end
endfunction

function [lambda,x1,k,n_erro]=Metodo_potencia_2(A,x0,epsilon,M)
    
k=1
x0=x0/norm(x0,2)    
x1=A*x0
while k<=M
    lambda=x0'*x1 //Quociente de Rayleigh, x0 e unitario
    if lambda<0
        x1=-x1 //mantém x1 com o mesmo sentido de x0
    end
    x1=x1/norm(x1,2)
    n_erro=norm(x1-x0,2)
    if n_erro < epsilon
        disp("O metodo convergiu com a precisao passada")
        break
    end    
    x0=x1
    x1=A*x0    
    k=k+1
end

if k>M    
    disp("O numero de iteracoes passou do maximo de iteracoes permitido")
end
endfunction
