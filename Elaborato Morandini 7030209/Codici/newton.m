function [x , it,count] = newton (f , f1 , x0 , maxIt , tol )
% Determina uno zero della funzione f applicando il metodo di Newton .
% Input :
% f - Funzione di cui voglio trovare la radice,  f1 - derivata prima della funzione f, x0 - approssimazione iniziale della radice
% maxIt - numero massimo di iterazioni [DEFAULT 100], tol - tolleranza [DEFAULT 10^-3]
% Output :  x - approssimazione della radice  it - numero di iterazioni fatte  count - numero di valutazioni funzionali fatte
if maxIt <0
     maxIt=100;
end
if tol<0
     tol=1e-3;
end
count=0;
x = x0 ;
for i =1: maxIt
     x0 = x ;
     fx = feval (f , x0 );
     f1x = feval ( f1 , x0 );
     count=count+2;
     if f1x==0 
         break
     end
     x = x0 - fx / f1x ;
     if abs (x - x0 ) <= tol *(1+ abs ( x ))
         break
     end
end
it = i ;
if ( abs (x - x0 ) > tol *(1+ abs ( x )))
     disp (" il metodo non converge ") 
end
end