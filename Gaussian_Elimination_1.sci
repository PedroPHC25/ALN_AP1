//////////////////////////////////////////////////////////////////////////
//Variáveis de saída:
//x: solução do sistema Ax=b (assumimos que tal solução existe).
//C: Seja A=LU a decomposição LU de A.
//Então C(i,j)=L(i,j) para i>j e C(i,j)=U(i,j) para j>=i.
//////////////////////////////////////////////////////////////////////////
function [x, C]=Gaussian_Elimination_1(A, b)
    C=[A,b];
    [n]=size(C,1);
    for j=1:(n-1)
        //O pivô está na posição (j,j)
        for i=(j+1):n
            //O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
            C(i,j)=C(i,j)/C(j,j);
            //Linha i  Linha i - C(i,j)*Linha j
            //Somente os elementos da diagonal ou acima dela são computados
            //(aqueles que compõem a matrix U)
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end
    x=zeros(n,1);
    // Calcula x, sendo Ux=C(1:n,n+1)
    x(n)=C(n,n+1)/C(n,n);
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i+1:n)*x(i+1:n))/C(i,i);
    end
    C=C(1:n,1:n);
endfunction

/////////////////////////// EXERCÍCIO 1 ////////////////////////////////

disp("-------------------- Exercício 1 --------------------")

A_teste_1 = [1 4; 8 5]
b_teste_1 = [26; 46]
x_teste_1 = Gaussian_Elimination_1(A_teste_1, b_teste_1)

disp("Teste 1: Resultado esperado: [2; 6]. Resultado obtido:", x_teste_1)

A_teste_2 = [7 5 1; 0 6 8; 9 3 4]
b_teste_2 = [45; 44; 40]
x_teste_2 = Gaussian_Elimination_1(A_teste_2, b_teste_2)

disp("Teste 2: Resultado esperado: [2; 6; 1]. Resultado obtido:", x_teste_2)

/////////////////////////// EXERCÍCIO 2 ////////////////////////////////

disp("-------------------- Exercício 2 --------------------")

A1 = [1 -2 5 0; 2 -4 1 3; -1 1 0 2; 0 3 3 1]
b1 = [1;0;0;0]
x1 = Gaussian_Elimination_1(A1, b1)

disp("x1:", x1)

/////////////////////////// EXERCÍCIO 3 ////////////////////////////////

disp("-------------------- Exercício 3 --------------------")

function [x, C]=Gaussian_Elimination_2(A, b)
    C=[A,b];
    [n]=size(C,1);
    for j=1:(n-1)
        // Se o pivô for 0, permuta as linhas
        if C(j,j) == 0 then
            // Encontrando a primeira linha sem zero abaixo do pivô
            row_without_zero = find(C(j+1:n,j),1)
            // Trocando essas linhas
            C([j j+row_without_zero],:) = C([j+row_without_zero j],:)
        end
        //O pivô está na posição (j,j)
        for i=(j+1):n
            //O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
            C(i,j)=C(i,j)/C(j,j);
            //Linha i  Linha i - C(i,j)*Linha j
            //Somente os elementos da diagonal ou acima dela são computados
            //(aqueles que compõem a matrix U)
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end
    x=zeros(n,1);
    // Calcula x, sendo Ux=C(1:n,n+1)
    x(n)=C(n,n+1)/C(n,n);
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i+1:n)*x(i+1:n))/C(i,i);
    end
    C=C(1:n,1:n);
endfunction

x1_2 = Gaussian_Elimination_2(A1, b1)

disp("x1_2:", x1_2)

A2 = [0 10^(-20) 1; 10^(-20) 1 1; 1 2 1]
b2 = [1; 0; 0]
x2 = Gaussian_Elimination_2(A2, b2)

disp("x2:", x2)













