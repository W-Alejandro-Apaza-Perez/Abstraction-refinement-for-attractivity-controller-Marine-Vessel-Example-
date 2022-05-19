function L=coordbox2state(Imin,Imax,p_x)

L=[0];
for k=length(Imin):-1:1
    Lp=[];
    for j=L'
        Lp=[Lp;(([Imin(k):Imax(k)]'-1)*p_x(k)+j)];
    end
    L=Lp;
end
L=L+1;