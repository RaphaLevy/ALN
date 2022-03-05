function [x,C]=Gaussian_Elimination_3(A,b)
  
C=[A,b];  
[n]=size(C,1);

for j=1:(n-1)
    //O pivô está na posição (j,j)
        //Se o pivô for igual a zero as linhas devem ser trocadas para chegar até o pivô máximo
        if C(j,j)==0 then
            [a, b] = max(abs(C(j:n, j)))
             C([j, b+(j-1)], :) = C([b+(j-1), j], :)
            end

        
        //O elemento C(i,j) é o elemento na posição (i,j) of L 
        //na decomposição LU de A
        for i=(j+1):n
            C(i,j)=C(i,j)/C(j,j);
        //Linha i <- Linha i - C(i,j)*Linha j
        //Somente os elementos da diagonal ou acima da diagonal são computados
        //(aqueles que compõem a matrix U)
           C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
    end
end

x=zeros(n,1);

// Calcula x, sendo Ux=C(1:n,n+1)

x(n)=C(n,n+1)/C(n,n);
for i=n-1:-1:1
    x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
end

C=C(1:n,1:n);

endfunction
