 %%Inverted Pendulum 
clear
syms x y a
b=1; %value of b
f=[y ; sin(x)-a*x*cos(x)+y*cos(x)]; % ODE equations
J=jacobian(f,[x y a]) % Find jacobian
eq=[0 ; 0 ; 1] %The equilibrium 
J0=subs(J,[x;y;a],eq) %insert the equilibrium in J
N=null(J0) % Null space
%%
syms y1 y2
alpha=(sin(y1)-((y2+1)*y1)*cos(y1))/(1+cos(y1))
t=taylor(alpha,[y1,y2],[0  0],'order',3) % taylor expanison 
M=subs(hessian(alpha),[y1,y2],[0  0]);% hessian is a square matrix of second-order partial derivatives 
% of a scalar-valued function, subs subtitute
[U,D]=eig(M) % eig(M) is eigenvalue and U is eigenverctor
u1=U(:,1) %first eigenverctor 
u2=U(:,2) % second eigenverctor
S=M-D(1,1)*eye(2) %  to find the null space of M, where eye(2)is identity matrix and D(1,1) is first eigenvalue (M-lambda I)
null(S)% nullspace of S, we get the eigenvector
B1=diff(alpha,y1) %first derivative w.r.t y1
B2=diff(alpha,y2) %first derivative w.r.t y2
C11=diff(B1,y1) %second derivative w.r.t y1
simplify(C11)
C12=diff(B1,y2) %second derivative w.r.t y2
simplify(C12)
C21=diff(B2,y1) %second derivative w.r.t y1
simplify(C21)
C22=diff(B2,y2) %second derivative w.r.t y2
simplify(C22)

