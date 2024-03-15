function yy = hermite(xi,fi, f1i, xx)
%Input
%xi: Vettore colonna contenete le ascisse d'interpolazione distinte, 
%fi: Vettore colonna contenete i valori della funzione nelle ascisse d'interpolazion
%f1i: Vettore colonna contenente i valori che la derivata prima della funzione assume nelle ascisse d'interpolazione
%xx: Vettore contenente le ascisse in cui vogliamo approssimare la funzione
%Output
%yy: Valori approssimati della funzione con il polinomio interpolante di Hermitte in forma di Newton
%Calcola i valori approssimati della funzione(di cui conosciamo sia i valori Y che assume nelle ascisse X ed i valori Y1  
%la cui derivata prima assume nelle stesse ascisse) calcolati attraverso il polinomio interpolante in
%forma di Lagrange nelle ascisse XQ.

if(isempty(xx)) 
    error("Il vettore contenente le ascisse in cui interpolare la funzione è vuoto!"),
end
if(length(xi)~=length(fi))
    error("Numero di ascisse d'interpolazione e di valori della funzione non uguale!"),
end
if (length(xi) ~= length(unique(xi)))
    error("Ascisse d'interpolazione non distinte!"),
end %uso la function unique che restuisce un vettore contenente i valori senza ripetizioni di X
if(length(f1i)~=length(fi))
    error("Lunghezza dei dati in ingresso non compatibile!"),
end
if(size(xi,2)>1||size(fi,2)>1||size(f1i,2)>1)
    error("Inserire vettori colonna!"),
end
n=(length(xi));
fi(1:2:2*n-1)=fi;
fi(2:2:2*n)=f1i;
df=diffdivher(xi,fi');
n=length(df)-1;
yy=df(n+1)*ones(size(xx));%algoritmo di horner per il calcolo dei valori di un polinomio
for i=n:-1:1
    yy=yy.*(xx-xi(ceil(i/2)))+df(i);
end
return
end

function f=diffdivher(x,f)
%function per calcolare le differenze divise per il polinomio interpolante di hermite
n=(length(f)/2)-1;
for i=2*n+1:-2:3
    f(i)=(f(i)-f(i-2))/(x(ceil(i/2))-x(ceil((i-1)/2)));
end
for j= 2:2*n+1
    for i=(2*n+2):-1:j+1
        f(i)=(f(i)-f(i-1))/(x(ceil(i/2))-x(ceil((i-j)/2)));
    end
end
end

