function x = mialdl(A,b)
%Calcola la soluzione del sistema lineare Ax=b utilizzando la fattorizzazione LDLT
% Input:: A = matrice sdp da fattorizzare, b = vettore termini noti
% Output: x = soluzione del sistema
[m,n] = size(A);
if m~=n
    error("Errore: La matrice deve essere quadrata");
end
if A(1,1)<=0
    error('Errore: La matrice deve essere sdp');
end
A(2:n,1)=A(2:n,1)/A(1,1); %fattorizzazione LDLT
for i=2:n
    v = (A(i,1:i-1).').*diag(A(1:i-1,1:i-1));
    A(i,i) = A(i,i)-A(i,1:i-1)*v;
    if A(i,i)<=0
        error("Errore: La matrice deve essere sdp");
    end
    A(i+1:n,i) = (A(i+1:n,i)-A(i+1:n,1:i-1)*v)/A(i,i);
end
x=b;
for i=1:n
    x(i+1:n) = x(i+1:n)-(A(i+1:n,i)*x(i));
end
x = x./diag(A);
for i=n:-1:2
    x(1:i-1) = x(1:i-1)-A(i,1:i-1).'*x(i);
end
end
