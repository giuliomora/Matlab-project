function ll = lebesgue(a, b, nn, type)
% Input: a,b = inizio e fine intervallo, nn = vettore con specificato il grado dei polinomi, type = specifica che tipo di ascisse di interpolazione usare
% se 0 utilizza le ascisse equispaziate nell'intervallo [a,b] altrimenti se 1 utilizza le ascisse di chebyshev
% Output: ll = vettore delle costanti di Lebesgue per ogni grado specificato in input
if a >= b
    error('Errore: a deve essere minore di b.');
end
if any(mod(nn, 1) ~= 0) || any(nn < 0)
    error('Errore: nn deve contenere solo numeri interi non negativi.');
end
if type ~= 0 && type ~= 1
    error('Errore: type deve essere 0 o 1.');
end
ll = ones(numel(nn), 1);
xq = linspace(a, b, 10001);
for j = 1:numel(nn)
    if type == 0 % Ascisse equidistanti
        x = linspace(a, b, nn(j)+1);
else % Ascisse di Chebyshev
    x = chebyshev(nn(j), a, b);
    end
lin = lebesgue_function(x,xq);
ll(j) = max(abs(lin));
end
end
function lin = lebesgue_function(x,xq)
n = length(x);
if length(unique(x)) ~= n
    error('Errore: Alcune ascisse sono uguali.');
end
lin = zeros(size(xq));
for j = 1:n
    L = ones(size(xq));
    for i = 1:n
        if i ~= j
            L = L .* (xq - x(i)) / (x(j) - x(i));
        end
    end
lin = lin + abs(L);
end
end

