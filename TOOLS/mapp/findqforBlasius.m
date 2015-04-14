format long;
options = odeset('AbsTol',1e-9,'RelTol',1e-6);

for i=1:10000 q(i)=0.33+0.01/10000*i;
[t,x]=ode45(@BlasiusFunc, [0,10],[0 0 q(i)]);
r(i)=x(length(x),2); 
% x(length(x),2) is the value of f2(10)
end;
plot(q,r);grid;xlabel('q');ylabel('r');