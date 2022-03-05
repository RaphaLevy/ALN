function [] = visualization(A) 
[m,n]=size(A);

subplot(1,1,1);imshow(A) // subplot(1,1,1) cria uma janela ("matriz") de tamanho 1x1, 
                         // e adiciona a figura na única posição disponível. A janela 
                         // é então mostrada com imshow(A)

endfunction
