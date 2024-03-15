function [I, nfeval] = newtonCotesAdatt(fun, a, b, tol, fa, f2, fm, f3, fb)
% Input: fun = funzione integranda, a = estremo inferiore dell'intervallo di integrazione, 
% b = estremo superiore dell'intervallo di integrazione
% tol = tolleranza richiesta, [fa,f2,f4,f5,fm] = valori calcolati negli estremi di integrazione e nei punti, intermedi. 
% Sono parametri opzionali che %diminuiscono il numero di valutazioni necessarie
% Output: I: approssimazione dell'integrale, nfeval = numero di valutazioni funzionali effettuate
if a > b
    error("Errore: gli estremi dell'intervallo non vanno bene");
end
if tol < 0
    error("Errore: la tolleranza specificata è minore o uguale a 0");
end
xm = (a+b)/2;
x2 = (a+xm)/2;
x3 = (xm+b)/2;
nfeval = 0;
if nargin == 4
    fa = feval(fun, a);
    fb = feval(fun, b);
    fm = feval(fun, xm);
    f2 = feval(fun, x2);
    f3 = feval(fun, x3);
    nfeval = nfeval + 5;
end
h = (b-a)/90;
x4 = (a+x2)/2;
x5 = (x2 + xm)/ 2;
x6 = (xm + x3)/ 2;
x7 = (x3 + b) / 2 ;
f4 = feval(fun, x4);
f5 = feval(fun, x5);
f6 = feval(fun, x6);
f7 = feval(fun, x7);
nfeval = nfeval + 4;
Itemp = h * (7*fa + 32*f2 + 12*fm + 32*f3 + 7*fb);
I = h/2 * (7*fa + 32*f4 + 12*f2 + 32*f5 + 14*fm + 32*f6 + 12*f3 + 32*f7 + 7*fb);
err = abs(I - Itemp) / 63;
if err > tol
    [I1, nfeval1] = NewtonCotesAdatt(fun, a, xm, tol / 2, fa, f4, f2, f5, fm);
    [I2, nfeval2] = NewtonCotesAdatt(fun, xm, b, tol / 2, fm, f6, f3, f7, fb);
    I = I1 + I2;
    nfeval = nfeval + nfeval1 + nfeval2;
end
end

