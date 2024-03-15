function rad=radice(x)
%rad=radice(x)
%Input: x, numero positivo di cui si vuole conoscere la radice
%Output: rad, radice quadrata del numero in input
%Restituisce la radice quadrata del numero x passato in input.
if (x==0)
    rad=0;
    return;
end
if x<0, error("x non puo' essere negativo"), end
x0=x;
for i=0:10000
    rad=x-((x^6 -x0)/(6*x^5));%passo iterativo
    if(abs(rad-x)<eps*(1+abs(x))) %controllo di aver raggiunto la massima precisione possibile
        break;
    end
    x=rad;
end
return;

