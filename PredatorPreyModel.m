%Example PredatorPreyModel in my second year report this example from Strogatz-Nonlinear-Dynamics and chaos
% This solution for 2 dim
clear
syms a b x y
eq1=a-x-(4*x*y)/(1+x^(2)); %First equation (ODE RHS)
eq2=b*x*(1-y/(1+x^(2))); %Second equation (ODE RHS)
assume([x>0,y>0,a>0,b>0])
J=jacobian([eq1,eq2],[x,y]); %find jacobian 
S1=solve([eq1,eq2],[x,y]); %find the values of x and y (Equilibrium)
J0=subs(J,S1) %insert the values of x and y in J 
bhopf=solve(trace(J0),b) % From trace we find the value of b (Hpof)
det(J0) % determinant of matrix J0
subs(det(J0),a,10) %replace all values of a by 10 in det(J0)
syms xh yh bh in my report delta_x,delta_y,delta_b (small)
H=subs([eq1;eq2],[x  y  b],[S1.x+xh,S1.y+yh,bhopf+bh]) %insert [S1.x+xh,S1.y+yh,bhopf+bh into x  y  b in 2 equations (Expanison)
H0=subs(H,a,10)%%replace all values of a by 10 
taylor(H0(1),[xh,yh,bh],'order',3)
taylor(H0(2),[xh,yh,bh],'order',3) 
%%
clear
syms x y a b s T
eq1=T*(a-x-(4*x*y)/(1+x^(2))); %First equation ODE RHS,x and y are functions a ,b are parameters
eq2=T*(b*x*(1-y/(1+x^(2))));
syms x0 y0 xc yc omega
k=1; % number of cos & sin components in harmonic expansion of x
xc=sym('xc',[1,k]);
yc=sym('yc',[1,k]);% unknown harmonic coefficients (cos)
xs=sym('xs',[1,k]); % unknown harmonic coefficients (sin)
xe=x0+cos((1:k)*s)*xc+sin((1:k)*s)*xs;% expansion of x(t)
ye=y0+cos((1:k)*s)*yc
dxe=diff(xe,s); % x'(t)=
dye=diff(ye,s);%
a=10
syms x1 y1 T1 b1 xh yh Th bh
vars=[x0;xc(1);y0;yc(1);xs(1); Th         ;bh] % list of unknowns
Ularge=[x;y;T;b] % list of unknowns, [x(t),y(t),b](large)
Usmall=[xh;yh;Th;bh] % list of unknowns (small)
Uhopf=  [2; 5  ; 1/sqrt(7) ;7/2] %Ularge-Usmall=Uhopf,Uhopf [x,y,T,b]=[a/5,a^2/25,1/sqrt(omega),3a/5-25/a] ,(delta(x),delta(y) are small),s0=[xhpof;0;yhpof;0;0;Thopf,bhopf] 
%SU=subs([eq1  eq2],Ularge,[2+x1  5+y1  T1+(1/sqrt(sym(7)))    (7/2)+b1 ]) %notice Uhopf+Usmall
SU=subs([eq1  eq2],Ularge,Uhopf+Usmall) %notice Uhopf+Usmall
TA=taylor(SU,Usmall,0*Usmall,'order',3) %
deq=subs(TA,[xh  yh],[xe,ye]) %[xe,ye]=[x0 + xc1*cos(s) + xs1*sin(s), y0 + yc1*cos(s)],where x0,xc1,xs1,yc1 are small(amplitude)
S=[1;reshape([cos(s*(1:k)),sin(s*(1:k))],[],1)]
S1=simplify(S*[deq(1)-dxe , deq(2)-dye]) %dxe,dye differentail eq
aeq=(int(S1,s,0,2*pi))/(2*pi)
S2=simplify(aeq)
TT=taylor(S2,[x0,xc(1),y0,yc(1),xs(1),T1,b1],[0 0 0 0 0 0 0],'order',2) %delta(t), TT is 6 equations
sgvals=[0;0;0;0;0;0;0] % valu
residual=subs(aeq,vars,sgvals) %
 J0=jacobian(TT(:),vars)
TT1=taylor(J0,[x0,xc(1),y0,yc(1),xs(1),Th,bh],[0 0 0 0 0 0 0],'order',2)
% J0=jacobian(aeq(:),vars)
J1=subs(J0,vars,sgvals)% Jaocbian is singular
V=null(J1)
V' %transpose
W=null(J1') 
%%
z=sym('z',size(vars))
yhat=sym('yhat',[size(V,2),1])%why we choose equations 
alpha=sym('alpha',[size(W,2),1])
EQ1=subs(S2(:),vars,V*yhat+z)+W*alpha %Applying equation 59 f(V*yhat+z)+W*alpha in equations after integral 
TEQ1=taylor(EQ1,[z;alpha;yhat],zeros(length(z)+length(alpha)+length(yhat),1),'order',2)
EQ2=V'*z %equation 60 V'*z=0
EQ=[EQ1;EQ2] %both equations
J3=jacobian(EQ,[z;alpha]) %jacobian both equations with respect to z & alpha
J01=subs(J3,[z;alpha;yhat],zeros(length(z)+length(alpha)+length(yhat),1)) %we can use [0*z;0*alpha;0*yhat]
% zeros(length(z)+length(alpha)+length(yhat),1)
R=EQ-J01*[z;alpha] %R is reminder 
TR=taylor(R,[z;alpha;yhat],zeros(length(z)+length(alpha)+length(yhat),1),'order',2) %R is small (1 is column)
S4=subs(-J01\R,[z;alpha],zeros(length(z)+length(alpha),1)) %first order iteration z1 alpha1 -S3*R^-1(inverse R),S4 is values of z and alpha
S5=subs(-J01\R,[z;alpha],S4)%insertaing S4 (second iteration )
TS=taylor(S5,[yhat],[0;0;0],'order',4)
sm=simplify(TS(8:9)/yhat(1)) %divid by yhat1
J02=jacobian(sm,[yhat(2),yhat(3)])
S6=subs(J02,[yhat(2);yhat(3)],[0 ; 0])
R=sm-S6*[yhat(2);yhat(3)]
Rsm=taylor(R,[yhat(2);yhat(3)],[0;0],'order',2) %R is small
S7=subs(-S6\R,[yhat(2);yhat(3)],[0;0])
%%



                                                        
