function x = mialu(A,b)
% x = mialu(A,b) 
% Data in ingresso una matrice A ed un vettore b, calcola la soluzione del sistema lineare Ax=b con il metodo di fattorizazione 
% LU con pivoting parziale 
% Input: A = matrice dei coefficienti, b = vettore dei termini noti 
% Output: x = soluzione del sistema lineare
[m,n] = size(A); 
if m ~= n 
    error("Errore: La matrice in input non è quadrata"); 
end 
if n~=length(b) 
    error("Errore: la dimensione del vettore b non coincide con la dimensione della matrice A"); 
end 
if size(b,2)>1 
error("Errore: il vettore b non è un vettore colonna"); 
end 
p = (1:n).'; 
for i=1:n 
    [mi,ki]=max(abs(A(i:n,i))); 
    if mi==0 
        error('Errore: La matrice non è non singolare.'); 
    end 
    ki = ki+i-1; 
    if ki>i 
        A([i,ki],:) = A([ki,i],:); 
        p([i,ki]) = p([ki,i]); 
    end 
    A(i+1:n,i) = A(i+1:n,i)/A(i,i); 
    A(i+1:n,i+1:n) = A(i+1:n,i+1:n)-A(i+1:n,i)*A(i,i+1:n); 
end 

x = b(p); 
for i=1:n 
    x(i+1:n) = x(i+1:n)-A(i+1:n,i)*x(i); 
end 
for i=n:-1:1 
    x(i) = x(i)/A(i,i); 
    x(1:i-1) = x(1:i-1)-A(1:i-1,i)*x(i); 
end 
end 
