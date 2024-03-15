function yy = myspline ( xi , fi , xx , type)
% Calcola la spline cubica interpolante naturale o not a knot a seconda del valore "tipo"
% Input: xi = vettore delle ascisse di interpolazione, fi = vettore dei valori della funzione nelle ascissedi interpolazione
% xq = punti in cui si vuole calcolare il polinomio interpolante, tipo = tipo di spline, 0 per naturale oppure 1 per not a knot
% Output: yy = valori che la spline assume nei punti xx
if nargin < 4 
    error ("Errore: numero degli argomenti errato") ;
elseif length ( xi ) ~= length ( fi ) 
    error ("Errore: i due vettori devono avere la stessa dimensione") ;
elseif length ( xi ) ~= length ( unique ( xi ) ) 
    error ("Errore: le ascisse devono essere distinte tra loro" ) ;
elseif size ( xi , 2) > 1 || size ( fi , 2) > 1 
    error ("Errore: vettore colonna non valido" )
elseif isempty ( xx ) 
    error ("Errore: partizione assegnata vuota") ;
end
n = length ( xi ) - 1;
epsilon = zeros (n -1 ,1) ;
q = zeros (n -1 ,1) ;
for i = 2 : n
    hi = ( xi ( i ) - xi (i -1) ) ;
    hi1 = ( xi ( i +1) - xi ( i ) ) ;
    q (i -1) = hi / ( hi + hi1 ) ;
    epsilon (i -1) = hi1 / ( hi + hi1 ) ;
end
l = size ( xi ) ;
diff_div = fi ;
l = l -1;
for j = 1 : 2
    for i = l +1 : -1 : j +1
        diff_div ( i ) = ( diff_div ( i ) - diff_div (i -1) ) /( xi ( i ) - xi (i - j ) ) ;
    end
end
diff_div = diff_div (3: l +1) ; % Calcolo le differenze divise
a = 2* ones (n -1 ,1) ;
if type == 0
    m = tridia (a , q , epsilon , diff_div * 6) ;
    m = [0; m ; 0];
else
    diff_div = diff_div * 6;
    a(1) = 2 - q(1);
    epsilon(1) = epsilon(1) - q(1);
    diff_div_1 = diff_div(1);
    diff_div(1) = (1 - q(1)) * diff_div(1);
    q(1) = 0;
    a(n-1) = 2 - epsilon(n-1);
    q(n-1) = q(n-1) - epsilon(n-1);
    diff_div_n = diff_div(n-1);
    diff_div(n-1) = (1 - epsilon(n-1)) * diff_div(n-1);
    epsilon(n-1) = 0;
    m = tridia(a, epsilon, q, diff_div);
    m0 = diff_div_1 - m(1) - m(2);
    mn = diff_div_n - m(n-1) - m(n-2);
    m = [m0; m; mn];
end
yy = zeros ( length ( xx ) , 1) ;
for j = 1 : length ( xx )
    for i = 2 : length ( xi )
        if (( xx ( j ) >= xi (i -1) && xx ( j ) <= xi ( i ) ) || xx ( j ) < xi (1) )
            hi = xi ( i ) - xi (i -1) ;
            ri = fi (i -1) - hi ^2/6* m (i -1) ;
            qi = ( fi ( i ) - fi (i -1) ) / hi - hi /6*( m ( i ) -m (i -1) ) ;
            yy ( j ) =(( xx ( j ) - xi (i -1) ) ^3* m ( i ) +( xi ( i ) - xx ( j ) ) ^3* m (i -1) ) /(6* hi ) + qi*( xx ( j ) - xi (i -1) ) + ri ;
            break
        end
    end
end
return ;
end

%Il codice utilizza la function di supporto tridia(a,b,c,g):

function x = tridia (a , b , c , g )
% Risolve il sistema lineare fattorizzabile LU
% Input :
% a = Vettore diagonale principale
% b = Vettore sovradiagonale
% c = Vettore sottodiagonale
% g = Vettore dei termini noti
% Output: x = Vettore soluzione
n = length ( g ) ;
x = g ;
for i = 1 : n -1
    b ( i ) = b ( i ) / a ( i ) ;
    a ( i +1) = a ( i +1) - b ( i ) * c ( i ) ;
    x ( i +1) = x ( i +1) - b ( i ) * x ( i ) ;
end
x ( n ) = x ( n ) / a ( n ) ;
for i = n -1 : -1 : 1
    x ( i ) = ( x ( i ) - c ( i ) * x ( i +1) ) / a ( i ) ;
end
return ;
end

