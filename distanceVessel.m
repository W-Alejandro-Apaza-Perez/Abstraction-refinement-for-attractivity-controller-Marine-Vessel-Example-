function K_out=distanceVessel(x,d,bound,xref,obstacle1,obstacle2)
tol=10^-6; %tolerance

x1_min=x(1)-d(1);
x1_max=x(1)+d(1);
x2_min=x(2)-d(2);
x2_max=x(2)+d(2);
x3_min=x(3)-d(3);
x3_max=x(3)+d(3);

if bound(1,1)<=(x1_min+tol) && x1_max<=(bound(1,2)+tol) && bound(2,1)<=(x2_min+tol) && (x2_max+tol)<=bound(2,2) && bound(3,1)<=(x3_min+tol) && (x3_max+tol)<=bound(3,2)...
    && ((x1_max-tol)<=obstacle1(1,1) || (x2_max-tol)<=obstacle1(2,1) || (x1_min+tol)>=obstacle1(1,2) || (x2_min+tol)>=obstacle1(2,2))...
    && ((x1_max-tol)<=obstacle2(1,1) || (x2_max-tol)<=obstacle2(2,1) || (x1_min+tol)>=obstacle2(1,2) || (x2_min+tol)>=obstacle2(2,2))
  K_out=norm(x(1:3)-xref(1:3)); 
else
  K_out=1000;
end