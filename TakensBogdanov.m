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
xc=sym('xc',[1,k]); % unknown harmonic coefficients (cos)
xs=sym('xs',[1,k]); % unknown harmonic coefficients (sin)
xe=x0+cos((1:k)*t)*xc.'+sin((1:k)*t)*xs.'; % expansion of x(t)
dxe=diff(xe,t); % x'(t)
d2xe=diff(dxe,t); %x''(t)
deq=subs(subs(omega^2*d2xe-f,[x,y],[xe,omega*dxe]),xs(1),0); %-om^2 x''+f(x,om x',a), fixing xs(1)=0
aeq=simplify(int([1;reshape([sin(t*(1:k));cos(t*(1:k))],[],1)]*deq,t,0,2*pi)/2/pi) % resulting system of algebraic eqs
vars=[x0;xc(1);xs(1);reshape([xs(2:end);xc(2:end)],[],1);omega;a] % list of unknowns
sgvals=[0;0;0;zeros(2*k-2,1);sqrt(2);0] % values of unknowns at singular point
residual=subs(aeq,vars,sgvals) % check if singular point is solution
J0=subs(jacobian(aeq,vars),vars,sgvals) % Jaocbian is singular
V=null(J0)
W=null(J0')
%%
syms y1 y2 y3 c x0 xc1 
A=subs(aeq,[x0,xc1,omega,a],[((1/sqrt(5)*y3)-(2*c)),y1,(sqrt(2)+y2),((2/sqrt(5)*y3)+c)])
%To find c,I will choose A(1)
A(1) %First equation
J1=jacobian(A(1),c)
SU=subs(J1,[c,y3],[0 0])
RHS=A(1)-SU*c %Rarrange the equation 
RL=-RHS/SU
RL1=simplify(subs(RL,c,RL))
T1=taylor(RL1,[y1,y3,c], [0  0  0],'order',3) % value of c
%Inserting T1(the value of c) into A(2) to get y3
A(2)%Second equation 
SU1=expand(subs(A(2),c,T1)/y1)
J11=jacobian(SU1,y3)
SU2=simplify(subs(J11,[y2,y3],[0  0]))
RHS1=simplify(SU1-SU2*y3)
RL2=simplify(-RHS1/SU2)
RL3=simplify(subs(RL2,y3,RL2))
T2=simplify(taylor(RL3,[y1,y2,y3],[0,0,0],'order',3)) %value ofy3
%Inserting T1(the value of c) into A(3)(second equation) to get y2
A(3) %Third equation
SU3=expand(subs(A(3),c,T1)/y1)
J2=jacobian(SU3,y2)
SU4=subs(J2,y2,0)
RHS2=simplify(SU3-SU4*y2)
RL4=-RHS2/SU4
RL5=simplify(subs(RL4,y2,RL4))
T3=simplify(taylor(RL5,[y1,y2],[0,0],'order',2)) 
RL6=subs(T3,y3,T2) 
T4=simplify(taylor(RL6,[y1,y2],[0,0],'order',3)) %the value of y2

