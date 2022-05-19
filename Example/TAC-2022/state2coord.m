function I=state2coord(i,p_x)

I=zeros(length(p_x)-1,1);
for k=1:length(p_x)-1
    I(k)=floor(mod(i-1,p_x(k+1))/p_x(k))+1; %floor rounds the number to -inf
end