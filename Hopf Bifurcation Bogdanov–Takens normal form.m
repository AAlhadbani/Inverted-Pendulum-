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
clear
syms x1 y1 b1
xdot=(7*x1/5)-(8*y1/5)+12*(((x1)*(y1))/25)-8*(x1)^2/25
ydot=(28*x1/5)-(7*y1/5)+8*(((b1)*(x1))/5)-2*(((b1)*(y1))/5)+21*(((x1)*(y1))/50)-7*(x1)^2/25
s1x=solve(xdot,x1)
s1y=solve(ydot,x1)
s2x=taylor(s1(1),[x1,y1],[0 0],'order',2)
s2y=taylor(s1(1),[x1,y1,b1],[0 0 0],'order',3)
s3=solve(s2,y1)
%%
clear
syms y2 y1 y3 
A =20*c^2 - 4*5^(1/2)*c*y3 - 25*c + y3^2 - (5*y1^2)/2
J1=jacobian(A,c)
%ans =40*c - 4*5^(1/2)*y3 - 25
J10=subs(J1,[c,y1,y3],[0 0 0])
RHS=A-J10*c
RL=-RHS/J10
RL1=simplify(subs(RL,c,RL))
T0=taylor(RL,[c,y1,y3],[0  0  0],'order',3)
T1=taylor(RL1,[c,y1,y3],[0  0  0],'order',3)
%%
A=((-sqrt(2)*y1^2)/10)+((sqrt(2)*y3^2)/25)-((y1^2*y2)/10)+((y2*y3^2)/25)-((y2*y3)/2*sqrt(5))-y3/sqrt(10)
 
A =(2^(1/2)*y3^2)/25 - (2^(1/2)*y1^2)/10 - (10^(1/2)*y3)/10 - (y1^2*y2)/10 + (y2*y3^2)/25 - (5^(1/2)*y2*y3)/2
M=solve(A,y3)
M =(5*(5*5^(1/2)*y2 + 10^(1/2) - 10*((2*y1^2*y2^2)/125 + (4*2^(1/2)*y1^2*y2)/125 + (4*y1^2)/125 + (5*y2^2)/4 + (5^(1/2)*10^(1/2)*y2)/10 + 1/10)^(1/2)))/(4*(y2 + 2^(1/2)))
(5*(5*5^(1/2)*y2 + 10^(1/2) + 10*((2*y1^2*y2^2)/125 + (4*2^(1/2)*y1^2*y2)/125 + (4*y1^2)/125 + (5*y2^2)/4 + (5^(1/2)*10^(1/2)*y2)/10 + 1/10)^(1/2)))/(4*(y2 + 2^(1/2)))
taylor(M(1),[y1,y2],[0 0],'order',3)
ans =
-(2^(1/2)*10^(1/2)*y1^2)/10
%%
    clear
syms y2 y1 y3 
A=((-sqrt(2)*y1^2)/10)+((sqrt(2)*y3^2)/25)-((y1^2*y2)/10)+((y2*y3^2)/25)-((y2*y3)/2*sqrt(5))-y3/sqrt(10)
J1=jacobian(A,y3)
J10=subs(J1,[y2,y3,y1],[0 0 0])
RHS=A-J10*y3
RL=-RHS/J10
RL1=simplify(subs(RL,y3,RL))
T0=taylor(RL,[y3,y1,y2],[0  0  0],'order',3)
T1=taylor(RL1,[y3,y1,y2],[0  0  0],'order',3)
\end{lstlisting}
%%
clear
syms y1 y2 y3
A =sym([0  -1/sqrt(10);-2^(1/2)   -1/sqrt(5)])%invertable matrix
AA=inv(A) %Inverse invertible matrix A
y2=-(sqrt(2)*y1^2)/5
y3=y1^2/sqrt(5)
B=[(-sqrt(2)*y1^2/10)+(sqrt(2)*y3^2/25)-(y1^2*y2/10)+(y2*y3^2/25)-(y2*y3/2*sqrt(5));(y1^2/5)+(2*y3^2/25)-(y2^2/2)] %orginal matrix
C=[- (2^(1/2)*y1^2)/10 - (4*2^(1/2)*y1^4)/125 - (2^(1/2)*y1^6)/1250;(13*y1^4)/500 + y1^2/5] %after inserting y2 and y3
%M=AA*B %multiply inverse matrix by B 
M=AA*C %multiply inverse matrix by C
S1=subs(M(1),[y2,y3],[sqrt(2)*y1^2/5  y1^2/sqrt(5)]) 
S2=subs(M(2),[y2,y3],[sqrt(2)*y1^2/5  y1^2/sqrt(5)])
TS1=taylor(S1,y1,0,'order',3) %find taylor S1, where S1 the first row in matrix M
TS2=taylor(S2,y1,0,'order',3) %find taylor S2, where S2 the second row in matrix M
%%
B=(y1^(2)/5)+(2*y3^(2)/25)+(y2^(2)/2)+(sqrt(2)*y2)+(y3/sqrt(5))
So1=solve(B,y3)
T=taylor(So1(1),[y1,y2],[0,0],'order',2)
 %%
clear
syms y2 y1 y3 

A=y1^2/5 + y2^2/2 + 2^(1/2)*y2 + (2*y3^2)/25 + (5^(1/2)*y3)/5;
J1=jacobian(A,y3)
J10=subs(J1,[y2,y3,y1],[0 0 0])
RHS=A-J10*y3
RL=-RHS/J10
RL1=simplify(subs(RL,y3,RL))
T0=taylor(RL,[y3,y1,y2],[0  0  0],'order',3)
T1=taylor(RL1,[y3,y1,y2],[0  0  0],'order',3)



