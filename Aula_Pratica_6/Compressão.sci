function [] = compression(A, p) 
[m,n]=size(A);    
if p < 0 | p > 1 //p é uma probabilidade, logo deve estar entre 0 e 100%
    break
end
AA = double(A)  //converte os elementos da matriz A para valores reais
[U,S,V]=svd(AA) //decomposição svd para uma matriz nativa do Scilab, 
                //S matriz diagonal de valores singulares de AA, 
                //U,V matrizes quadradas ortogonais ou unitárias de vetores singulares
r = size(S)(1)     //número de valores singulares positivos
s=max(1,int(p*r)) //pega os s maiores valores singulares de AA (s<r)
C=U(:,1:s)*S(1:s,1:s)*V'(1:s,:) //gera a imagem original e a comprimida lado a lado
C=iconvert(C,11) //converte para o formato inteiro com um byte de tamanho entre 0 e 225
subplot(1,2,1);imshow(A);subplot(1,2,2);imshow(C)

endfunction
