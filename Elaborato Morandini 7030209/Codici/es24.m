data = load('data.mat').data;
x = data(:,1);
y = data(:,2);
% Inizializza vettore errore
error = zeros(1,15);
for m = 1:15 %Calcola i coefficienti del polinomio di approssimazione ai minimi quadrati che meglio approssima i dati ordine decrescente di grado.
    p = polyfit(x, y, m); % Calcola il vettore di valori del polinomio approssimante sui punti x dei dati.
    y_fit = horner(p, x); % Calcola l'errore di approssimazione tra i dati originali e quelli approssimati
    error(m) = norm(x.^m - y_fit);
end

% semilogy
figure;
semilogy(1:15, error, 'b*-');
xlabel('Grado del polinomio (n)');
ylabel('Error');
title("Errore dell'approssimazione");
grid on;
function y = horner(p, x)
y = p(1);
for i = 2:length(p)
    y = y .* x + p(i);
end
end

