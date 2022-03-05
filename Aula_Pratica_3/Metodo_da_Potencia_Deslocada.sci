function [lambda,x1,k,n_erro]=Potencia_deslocada_inversa(A,x0,epsilon,alfa,M)
    
//VARIAVEIS DE ENTRADA    
//A: matriz nxn diagonalizavel
//x0: vetor NAO NULO a ser usado como aproximacao inicial do autovetor dominante    
//epsilon: precisao a ser usada no criterio de parada
//alfa: valor do qual se deseja achar o autovalor de A mais proximo
//M: numero maximo de iteracoes

//VARIAVEIS DE SAIDA
//lambda: autovalor de A mais proximo de alfa
//x1: autovetor unitario (norma 2) correspondente a lambda
//k: numero de iteracoes necessarias para a convergencia
//n_erro: norma 2 do erro

k=1
x0=x0/norm(x0,2)  //coord maior modulo x0  
n = length(x0)
while k<=M
    x1 = Gaussian_Elimination_4_AP3(A-alfa*eye(n,n),x0) // resolve (A-alfa*I)*x1=x0
    x1=x1/norm(x1,2)
    lambda=x1'*A*x1 //Quociente de Rayleigh, x1 e unitario
    if lambda<0
        x1=-x1 //mantem x1 com o mesmo sentido de x0
    end    
    n_erro=norm(x1-x0,2)
    if n_erro < epsilon
        disp("O metodo convergiu com a precisao passada") 
        break
    end    
    x0=x1
    k=k+1
end

if k>M    
    disp("O numero de iteracoes passou do maximo de iteracoes permitido")   
end
endfunction

function [lambda,x1,k,n_erro]=Potencia_deslocada_Rayleigh(A,x0,epsilon,alfa,M)

k=1
lambda=alfa
x0=x0/norm(x0,2)  //coord maior modulo x0  
n = length(x0)
while k<=M
    x1 = Gaussian_Elimination_4_AP3(A-lambda*eye(n,n),x0) // resolve (A-lambda*I)*x1=x0
    x1=x1/norm(x1,2)
    lambda=x1'*A*x1 //Quociente de Rayleigh, x1 e unitario  
    n_erro=norm(abs(x1)-abs(x0), 2) //Pega o valor absoluto dos vetores
    if n_erro < epsilon
        disp("O metodo convergiu com a precisao passada") 
        break
    end    
    x0=x1
    k=k+1
end

if k>M     
    disp("O numero de iteracoes passou do maximo de iteracoes permitido")
end
endfunction
