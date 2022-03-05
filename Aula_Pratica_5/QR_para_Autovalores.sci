function [S] = espectro(A, tol) //A deve ser simétrica e quadrada
    
[Q_k,R_k]=qr_GS(A) //o arquivo Metodo_Gram_Schmidt precisa ser executado para o funcionamento deste
A_k = R_k * Q_k    //A_k é triangular superior

while norm(diag(A_k) - diag(A),%inf) > tol
    A = A_k
    [Q_k,R_k]=qr_GS(A)
    A_k = R_k * Q_k 
    norm(diag(A_k) - diag(A),%inf)
end 

S = diag(A_k);

endfunction
