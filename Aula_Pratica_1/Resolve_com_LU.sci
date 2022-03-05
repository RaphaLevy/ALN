function [X]=Resolve_com_LU(C,B,P)
X=B 





function [y,C1]=Gaussian_Elimination_1(L,b)
  
C1=[L,b];  
[n]=size(C1,1);

L = eye(n,m) + tril(C, -1)
U = triu(C)

for j=1:(n-1)
    for i=(j+1):n
        C1(i,j)=C1(i,j)/C1(j,j);
        C1(i,j+1:n+1)=C1(i,j+1:n+1)-C1(i,j)*C1(j,j+1:n+1);
    end
  end
  
y=zeros(n,1);

// Calcula y, sendo Lx=b(1:n,n+1)

y(n)=C1(n,n+1)/C1(n,n);
for i=n-1:-1:1
    y(i)=(C1(i,n+1)-C1(i,i:n)*y(i:n))/C1(i,i);
end

C1=C1(1:n,1:n);

endfunction


function [x,C2]=Gaussian_Elimination_1(U,y)
  
C2=[U,y];  
[m]=size(C2,1);

L = eye(n,m) + tril(C, -1)
U = triu(C)

for j=1:(m-1)
    for i=(j+1):m
        C2(i,j)=C2(i,j)/C2(j,j);
        C2(i,j+1:m+1)=C2(i,j+1:m+1)-C2(i,j)*C2(j,j+1:m+1);
    end
  end
  
x=zeros(m,1);

// Calcula x, sendo Ux=y(1:n,n+1)

x(m)=C2(m,m+1)/C2(m,m);
for i=m-1:-1:1
    x(i)=(C2(i,m+1)-C2(i,i:m)*x(i:m))/C2(i,i);
end

C2=C2(1:m,1:m);

endfunction




endfunction
