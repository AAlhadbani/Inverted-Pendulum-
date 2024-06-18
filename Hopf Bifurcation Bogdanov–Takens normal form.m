%% LSR for Hopf along line b=1 for Takens-Bogdanov normal form in algebraic form
% see Kuznetsov, Elements of Applied Bifurcation Theory'98, p321, Thm 8.4
% for origin of Equation
clear
syms x y a b sg omega t
sg=-1; % set +1 or -1
b=1;
f=a-2*b*x+x^2+sg*x*y; % example f
syms x0
k=1; % number of cos & sin components in harmonic expansion of x
xc=sym('xc',[1,k])% unknown harmonic coefficients (cos)
xs=sym('xs',[1,k]) % unknown harmonic coefficients (sin)
xe=x0+cos((1:k)*t)*xc.'+sin((1:k)*t)*xs.'% expansion of x(t)
dxe=diff(xe,t) % x'(t)
d2xe=diff(dxe,t) %x''(t)
deq=subs(subs(omega^2*d2xe-f,[x,y],[xe,omega*dxe]),xs(1),0) %-om^2 x''+f(x,om x',a), fixing xs(1)=0
aeq=simplify(int([1;reshape([sin(t*(1:k));cos(t*(1:k))],[],1)]*deq,t,0,2*pi)/2*pi) % resulting system of algebraic eqs
vars=[x0;xc(1);reshape([xs(2:end);xc(2:end)],[],1);omega;a]% list of unknowns
sgvals=[0;0;zeros(2*k-2,1);sqrt(2);0];% values of unknowns at singular point
residual=subs(aeq,vars,sgvals);% check if singular point is solution
J0=subs(jacobian(aeq,vars),vars,sgvals)% Jaocbian is singular
V=null(J0)
W=null(J0') 
%%
syms y1 y2 y3 c x0 xc1 
A=subs(aeq,[x0,xc1,omega,a],[((1/sqrt(5)*y3)-(2*c)),y1,(sqrt(2)+y2),((2/sqrt(5)*y3)+c)])
%To find c,I will choose A(1)
A(1)
J1=jacobian(A(1),c)
SU=subs(J1,[c,y3],[0 0])
RHS=A(1)-J1*c %Rarrange the equation 
RL=-RHS/J1
RL1=simplify(subs(RL,c,RL))
T1=taylor(RL1,[y1,c,y3],[0  0  0],'order',3) % value of c
%Inserting T1(the value of c) into A(2) to get y3
SU1=subs(A(2),c,T1)/y1
SS1=solve(SU1,y3)
T2=taylor(SS1(2),[y1,y2,y3],[0,0,0],'order',3) %value ofy3
%Inserting T1(the value of c) into A(3)(second equation) to get y2

%%
syms y1 y2 y3
A =sym([0  1/sqrt(10);2^(1/2)   1/sqrt(5)])%invertable matrix
AA=inv(A)%Inverse invertable matrix
y2=(sqrt(2)*y1^2)/5
y3=-y1^2/sqrt(5)
B=[(-sqrt(2)*y1^2/10)-(sqrt(2)*y3^2/25)-(y1^2*y2/10)-(y2*y3^2/25)-(y2*y3/2*sqrt(5));(y1^2/5)+(2*y3^2/25)-(y2^2/2)] %orginal matrix
C=[- (2^(1/2)*y1^2)/10 - (4*2^(1/2)*y1^4)/125 - (2^(1/2)*y1^6)/1250;(13*y1^4)/500 + y1^2/5] %after inserting y2 and y3
%M1=AA*B
M2=AA*C
S1=subs(M(1),[y2,y3],[sqrt(2)*y1^2/5  -y1^2/sqrt(5)])
S2=subs(M(2),[y2,y3],[sqrt(2)*y1^2/5  -y1^2/sqrt(5)])
TS1=taylor(S1,y1,0,'order',3) %find taylor S1, where S1 the first row in matrix M
TS2=taylor(S2,y1,0,'order',3) %find taylor S2, where S2 the second row in matrix M
 %%



