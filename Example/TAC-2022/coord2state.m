function i=coord2state(I,p_x)

i=(I-1)'*p_x(1:length(I))+1;