function [x,it,count] = secanti(f,x0,x1,maxIt,tolx)
% [x,it] = secanti(f,x0,x1,maxIt,tolx)
% Calcola uno zero della funzione f usando il metodo delle secanti.
% Input: f - funzione di cui voglio determinare la radice, x0 - prima approssimazione iniziale della radice  x1 - seconda approssimazione iniziale della % radice, maxIt - numero massimo d'iterazioni [DEFAULT 100], tolx - tolleranza [DEFAULT 10^-3]
% Output: x - approssimazione della radice, it - numero di iterazioni eseguite, count - numero di valutazioni funzionali fatte
if maxIt <0
     maxIt=100; end
if tolx<0
     tolx=1e-3; end
count=0;
f0=feval(f,x0);
f1=feval(f,x1);
count=count+2;
if f1==0
    x=x1; return;
end
x=(x0*f1-x1*f0)/(f1-f0);
for i=1:maxIt
    if abs(x-x1)<= tolx*(1+abs(x1)) 
        break
    end
    x0=x1;
    f0=f1;
    x1=x;
    f1= feval(f,x);
    count=count+1;
    if f1==0 
        break
    end
    x=(x0*f1-x1*f0)/(f1-f0);
end
    it=i;
    if abs(x-x1)>tolx*(1+abs(x1))
           disp("Il metodo non converge");
    end

end