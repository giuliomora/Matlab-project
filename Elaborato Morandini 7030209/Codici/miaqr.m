function [x,nr] = miaqr(A,b)
%
% [x,nr] = miaqr(A,b)
%
% La funzione calcola la fattorizzazione QR del sistema lineare
% sovraddimensionato e restituisce, oltre alla fattorizzazione, la norma
% euclideiana del vettore residuo
%
% Input:
% A = matrice da fattorizzare
% b = vettore dei termini noti
%
% Output:
% x = soluzione del sistema Ax=b
% nr = norma del vettore residuo
%
[m,n] = size(A);
if(n>m)
    error('Errore: il sistema in input non è sovradeterminato.');
end
k = length(b);
if k~=m
    error('Errore: La dimensione della matrice e del vettore non coincidono.');
end
for i = 1:n
    a = norm(A(i:m,i),2);
    if a==0
        error('Errore: La matrice non ha rango massimo.');
    end
    if A(i,i)>=0
        a = -a;
    end
    v1 = A(i,i)-a;
    A(i,i) = a;
    A(i+1:m,i) = A(i+1:m,i)/v1;
    beta = -v1/a ;
    A(i:m,i+1:n) = A(i:m,i+1:n)-(beta*[1;A(i+1:m,i)])*...
        ([1;A(i+1:m,i)]'*A(i:m,i+1:n));
end
for i=1:n
    v = [1;A(i+1:m,i)];
    beta = 2/(v'*v);
    b(i:k) = b(i:k)-(beta*(v'*b(i:k)))*v;
end
for j=n:-1:1
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1)-A(1:j-1,j)*b(j);
end
x = b(1:n);
nr = norm(b(n+1:m));
end