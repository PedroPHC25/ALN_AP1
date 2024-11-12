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
            //(aqueles que compõem a matrix U)+
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

// Teste 1
A_teste_1 = [1 4; 8 5]
b_teste_1 = [26; 46]
x_teste_1 = Gaussian_Elimination_1(A_teste_1, b_teste_1)

disp("--- Teste 1 ---")
disp("A:", A_teste_1)
disp("x:", x_teste_1)
disp("b:", b_teste_1)
disp("Ax:", A_teste_1*x_teste_1)

// Teste 2
A_teste_2 = [7 5 1; 0 6 8; 9 3 4]
b_teste_2 = [45; 44; 40]
x_teste_2 = Gaussian_Elimination_1(A_teste_2, b_teste_2)

disp("--- Teste 2 ---")
disp("A:", A_teste_2)
disp("x:", x_teste_2)
disp("b:", b_teste_2)
disp("Ax:", A_teste_2*x_teste_2)

// Teste 3
A_teste_3 = [2 8 9 3; 6 3 1 1; 5 1 6 4; 7 5 7 3]
b_teste_3 = [95; 46; 54; 85]
x_teste_3 = Gaussian_Elimination_1(A_teste_3, b_teste_3)

disp("--- Teste 3 ---")
disp("A:", A_teste_3)
disp("x:", x_teste_3)
disp("b:", b_teste_3)
disp("Ax:", A_teste_3*x_teste_3)

/////////////////////////// EXERCÍCIO 2 ////////////////////////////////

disp("-------------------- Exercício 2 --------------------")

A1 = [1 -2 5 0; 2 -4 1 3; -1 1 0 2; 0 3 3 1]
b1 = [1; 0; 0; 0]
x1 = Gaussian_Elimination_1(A1, b1)

disp("x1:", x1)
disp("A1*x1:", A1*x1)
disp("b1:", b1)

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

// Testando de novo com A1 e b1
x1_2 = Gaussian_Elimination_2(A1, b1)

disp("Testando de novo com A1 e b1")
disp("x1_2:", x1_2)
disp("A1*x1_2:", A1*x1_2)
disp("b1:", b1)

// Testando com A2 e b2
A2 = [0 10^(-20) 1; 10^(-20) 1 1; 1 2 1]
b2 = [1; 0; 0]
x2 = Gaussian_Elimination_2(A2, b2)

disp("Testando com A2 e b2")
disp("x2:", x2)
disp("A2*x2:", A2*x2)
disp("b2:", b2)

/////////////////////////// EXERCÍCIO 4 ////////////////////////////////

disp("-------------------- Exercício 4 --------------------")

function [x, C]=Gaussian_Elimination_3(A, b)
    C=[A,b];
    [n]=size(C,1);
    for j=1:(n-1)
        // Se o pivô for 0, permuta com a linha com o maior pivô em módulo
        if C(j,j) == 0 then
            // Encontrando a linha com o maior pivô em módulo
            moduled_vector = abs(C(j+1:n,j))
            max_pivot_index = find(moduled_vector == max(moduled_vector))
            // Trocando essas linhas
            C([j j+max_pivot_index],:) = C([j+max_pivot_index j],:)
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

// Testando novamente com A2 e b2
x2_2 = Gaussian_Elimination_3(A2, b2)

disp("Testando novamente com A2 e b2")
disp("x2_2:", x2_2)
disp("A2*x2_2:", A2*x2_2)
disp("b2:", b2)

// Testando com A3 e b3
A3 = [10^(-20) 10^(-20) 1; 10^(-20) 1 1; 1 2 1]
b3 = b2

x3 = Gaussian_Elimination_3(A3, b3)

disp("Testando com A3 e b3")
disp("x3:", x3)
disp("A3*x3:", A3*x3)
disp("b3:", b3)

/////////////////////////// EXERCÍCIO 5 ////////////////////////////////

disp("-------------------- Exercício 5 --------------------")

function [x, C, P]=Gaussian_Elimination_4(A, b)
    C=[A,b];
    [n]=size(C,1);
    // Inicializando a matriz P
    P = eye(n,n)
    for j=1:(n-1)
        // Se o pivô não for o maior valor de sua coluna, troca as linhas
        if max(abs(C(j:n,j))) ~= abs(C(j,j)) then
            // Encontrando a linha com o maior pivô em módulo
            moduled_vector = abs(C(j:n,j))
            max_pivot_index = find(moduled_vector == max(moduled_vector))
            // Trocando essas linhas na C e na P
            C([j j+max_pivot_index-1],:) = C([j+max_pivot_index-1 j],:)
            P([j j+max_pivot_index-1],:) = P([j+max_pivot_index-1 j],:)
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

// Testando novamente com A3 e b3
x3_2 = Gaussian_Elimination_4(A3, b3)

disp("Testando novamente com A3 e b3")
disp("x3_2:", x3_2)
disp("A3*x3_2:", A3*x3_2)
disp("b3:", b3)

/////////////////////////// EXERCÍCIO 6 ////////////////////////////////

disp("-------------------- Exercício 6 --------------------")

function [X]=Resolve_com_LU(C, B, P)
    n = size(B,1)
    m = size(B,2)
    // Obtendo a L e a U da C
    L = eye(n,n) + tril(C,-1)
    U = triu(C)
    // Multiplicando dos dois lados por P para poder utilizar a decomposição LU de PA
    D = P*B
    // Resolvendo LY = D, onde Y = UX
    Y = zeros(n,m)
    Y(1,:) = D(1,:)
    for row = 2:n
        sub_factor = L(row,1:(row-1))*Y(1:(row-1),:)
        Y(row,:) = D(row,:) - sub_factor
    end
    // Resolvendo UX = Y
    X = zeros(n,m)
    X(n,:) = Y(n,:)/U(n,n)
    for row = (n-1):-1:1
        sub_factor = U(row,(row+1):m)*X((row+1):n,:)
        X(row,:) = (Y(row,:) - sub_factor)/U(row,row)
    end
endfunction

// Testando com A1 e B1
[x,C1,P1] = Gaussian_Elimination_4(A1,b1)

B1 = [2 4 -1 5; 0 1 0 3; 2 2 -1 1; 0 1 1 5]
X1 = Resolve_com_LU(C1,B1,P1)

disp("Testando com A1 e B1")
disp("X1", X1, "A1*X1", A1*X1, "B1", B1)

// Testando com A2 e B2
[x,C2,P2] = Gaussian_Elimination_4(A2,b2)

B2 = [1 1 2; 1 -1 0; 1 0 1]
X2 = Resolve_com_LU(C2,B2,P2)

disp("Testando com A2 e B2")
disp("X2", X2, "A2*X2", A2*X2, "B2", B2)

// Testando com A3 e B3
[x,C3,P3] = Gaussian_Elimination_4(A3,b3)

B3 = B2
X3 = Resolve_com_LU(C3,B3,P3)

disp("Testando com A3 e B3")
disp("X3", X3, "A3*X3", A3*X3, "B3", B3)
