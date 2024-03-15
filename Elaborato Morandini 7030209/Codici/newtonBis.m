function [x,nit] = newtonBis(fun,jacobian,x0,tol,maxit)
% Questo metodo risolve sistemi non lineare di equazioni attraverso l'uso del metodo di Newton. Può restituire inoltre il numero di iterazioni eseguite
% Input: fun = vettore delle funzioni, jacobian = matrice jacobiana di fun,  x0 = vettore delle approssimazioni iniziali, tol = tolleranza specificata
% maxit = quante iterazioni il metodo potrà fare al massimo, se non specificata maxit = 1000
% Output: x = vettore delle soluzioni, nit = numero di iterazioni svolte
if(nargin == 3)
    tol = eps;
end
if(nargin == 4)
    maxit = 1000;
end
if (tol<0)
    error('Errore: La tolleranza in input non può essere minore di 0.');
end
if (maxit<=0)
    error('Errore: Il massimo numero di iterazioni non può essere minore di 0.');
end
x = x0;
for nit=1:maxit
    x0 = x;
    fx0 = feval(fun,x0);
    fj = feval(jacobian,x0);
    x = x0-(fj\fx0);
    if (norm(x-x0)<=tol*(1+norm(x0)))
        disp('Tolleranza desiderata raggiunta.');
        break
    end
end
if not(norm(x-x0)<=tol*(1+norm(x0)))
    disp('Il metodo non è convergente.');
end
end