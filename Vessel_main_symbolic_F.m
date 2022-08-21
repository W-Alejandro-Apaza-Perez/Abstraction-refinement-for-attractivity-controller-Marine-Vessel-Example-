% Quantitative Reachability synthesis
clear
close all
clc
tol=10^-6; %tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Bounds
bound=[-3.5 3.5;-2.6 2.6;-pi pi];
inputbound=[-0.7 2;-1 1;-1 1];


obstacleState_1=[0.7 1;0.7 2.6;-pi pi];
obstacleState_2=[-1.8 -1.5;-2.6 -1;-pi pi];


% Target set

X_ref=[-0.5;-0.5;0];
X_refSet=[-1.2 0.1;-1.2 0.1;-0.5 0.5];
% Dynamics
%tol=10^-6;
T0 =1;

F = @(x, u) [x(1)+T0*(u(1)*cos(x(3))-u(2)*sin(x(3)));...
             x(2)+T0*(u(1)*sin(x(3))+u(2)*cos(x(3)));...
             x(3)+T0*u(3)];
         
A = @(u) [1 0 T0*(abs(u(1))+abs(u(2)));...
          0 1 T0*(abs(u(1))+abs(u(2)));...
          0 0 1];
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abstraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% Space discretization
n_x=[15;11;10];
p_x=[1;cumprod(n_x)];
n_states=p_x(end);                % total number of states
d_x=(bound(:,2)-bound(:,1))./n_x;    

% Input discretization
n_modes=[6;5;5];       % number of modes
d_u=(inputbound(:,2)-inputbound(:,1))./(n_modes-1);
mode=[];
for i=1:n_modes(1)
    for j=1:n_modes(2)
        for k=1:n_modes(3)
            mode=[mode inputbound(:,1)+[i-1;j-1;k-1].*d_u];
        end    
    end
end
n_modes=length(mode);

disp('Step 1: Computation of the abstraction''s transition relation')
h = waitbar(0,'Abstraction in progress...');

H=zeros(n_states,1); % H contains distance to safe set
Delta=zeros(n_states*n_modes,2);  % Delta contains the transition relation

for q=1:n_states 
    waitbar(q / n_states,h)
    c_q=state2coord(q,p_x); % coordinates of the cell
    x_c=bound(:,1)+(c_q-0.5).*d_x; % center of rectangle associated to state q %2x1 (in the actual state space)
    % Compute the distance to the safe set
    H(q)=distanceVessel(x_c,0.5*d_x,bound,X_ref,obstacleState_1,obstacleState_2);
    % Compute the transition relation
    for p=1:n_modes
        x_succ_c=F(x_c, mode(:,p));
        d_x_succ=0.5*A(mode(:,p))*d_x;
        x_succ=[x_succ_c-d_x_succ,x_succ_c+d_x_succ]; % Over-approximation of the reachable set from x under mode p
        if x_succ(3,1)-tol<-pi 
            x_succ(3,1)=x_succ(3,1)+2*pi;
        end
        if x_succ(3,2)-tol<-pi       %to cover all posibilities 
            x_succ(3,2)=x_succ(3,2)+2*pi;
        end
        if x_succ(3,2)+tol>pi 
            x_succ(3,2)=x_succ(3,2)-2*pi;
        end
        if x_succ(3,1)+tol>pi        %to cover all posibilities
            x_succ(3,1)=x_succ(3,1)-2*pi;
        end
        
         if bound(1,1)<=(x_succ(1,1)+tol) && (x_succ(1,2)-tol)<=bound(1,2) && bound(2,1)<=(x_succ(2,1)+tol) && (x_succ(2,2)-tol)<=bound(2,2)...
            && ((x_succ(1,2)-tol)<=obstacleState_1(1,1) || (x_succ(2,2)-tol)<=obstacleState_1(2,1) || (x_succ(1,1)+tol)>=obstacleState_1(1,2) || (x_succ(2,1)+tol)>=obstacleState_1(2,2))...
            && ((x_c(1)+0.5*d_x(1)-tol)<=obstacleState_1(1,1) || (x_c(2)+0.5*d_x(2)-tol)<=obstacleState_1(2,1) || (x_c(1)-0.5*d_x(1)+tol)>=obstacleState_1(1,2) || (x_c(2)-0.5*d_x(1)+tol)>=obstacleState_1(2,2))...
            && ((x_succ(1,2)-tol)<=obstacleState_2(1,1) || (x_succ(2,2)-tol)<=obstacleState_2(2,1) || (x_succ(1,1)+tol)>=obstacleState_2(1,2) || (x_succ(2,1)+tol)>=obstacleState_2(2,2))...
            && ((x_c(1)+0.5*d_x(1)-tol)<=obstacleState_2(1,1) || (x_c(2)+0.5*d_x(2)-tol)<=obstacleState_2(2,1) || (x_c(1)-0.5*d_x(1)+tol)>=obstacleState_2(1,2) || (x_c(2)-0.5*d_x(1)+tol)>=obstacleState_2(2,2))
            min_succ=coord2state(floor((x_succ(:,1)-bound(:,1))./d_x)+1,p_x); 
            max_succ=coord2state(ceil((x_succ(:,2)-bound(:,1))./d_x),p_x);
            Delta((q-1)*n_modes+p,:)=[min_succ,max_succ];   
        end
    end
end
close(h) 
toc

NC=[];
for q=1:n_states
         if max(~any(Delta((q-1)*n_modes+1:(q-1)*n_modes+n_modes,1)))
             NC=[NC;q];
         end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantitative synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
Hmax=max(H);%Inf; % Maximal distance to safe set

% Safety synthesis
Vn=H;
Vp=zeros(n_states,1);
C_mF=zeros(n_states,n_modes);
C=zeros(n_states,1);
Vmax=zeros(1,n_modes);
iter=0;

disp('Step 2: Safety Controller synthesis')
disp('Fixed point algorithm')
fixpoint=0;
while ~fixpoint
    iter=iter+1
    Vp=Vn;
    for q=1:n_states
        Vmax=zeros(1,n_modes);
        for p=1:n_modes
            if any(Delta((q-1)*n_modes+p,:)) % then we have successors
                Imin=state2coord(Delta((q-1)*n_modes+p,1),p_x);
                Imax=state2coord(Delta((q-1)*n_modes+p,2),p_x);
                if Imin(3)<=Imax(3)
                    Q_succ=coordbox2state(Imin,Imax,p_x)';             
                else
                    Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_x)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_x)'];
                end
                VN=[]; % This was the bug! need to clear this variable, because of the different lengths!
                for i=1:length(Q_succ)
                    VN(i)=Vp(Q_succ(i));
                end
                Vmax(p)=max(VN);  % Max value at successor states  
             else
                Vmax(p)=2*Hmax;
            end
        end
        [v,c]=min(Vmax);
        Vn(q)=max(H(q),v);  % MinMax value over inputs %NEW
        if any(Delta((q-1)*n_modes+c,:))
          C(q)=c;
        else
          C(q)=0;
        end
    end             
    fixpoint=isequal(Vp,Vn);
end

%Finding the controller for safety
for q=1:n_states
        for p=1:n_modes
            if any(Delta((q-1)*n_modes+p,:)) 
                Imin=state2coord(Delta((q-1)*n_modes+p,1),p_x);
                Imax=state2coord(Delta((q-1)*n_modes+p,2),p_x);
                if Imin(3)<=Imax(3)
                    Q_succ=coordbox2state(Imin,Imax,p_x)';             
                else
                    Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_x)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_x)'];
                end  
                VN=[]; 
                for i=1:length(Q_succ)
                    VN(i)=Vn(Q_succ(i)); 
                end
                if max(VN)<=Vn(q)  
                   C_mF(q,p)=1;  
                end    
            end     
        end
end


disp('Fixed point reached')
Iterations=iter

% Reachability synthesis
Wn=Vn;
Wp=zeros(n_states,1);
K=zeros(n_states,1);

iter=0;

disp('Step 3: Reachability Controller synthesis')
disp('Fixed point algorithm')
fixpoint=0;
while ~fixpoint
    iter=iter+1
    Wp=Wn;
    for q=1:n_states
        for p=1:n_modes
            if any(Delta((q-1)*n_modes+p,:)) 
                Imin=state2coord(Delta((q-1)*n_modes+p,1),p_x);
                Imax=state2coord(Delta((q-1)*n_modes+p,2),p_x);
                if Imin(3)<=Imax(3)
                    Q_succ=coordbox2state(Imin,Imax,p_x)';             
                else
                    Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_x)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_x)'];
                end
                VN=[]; 
                for i=1:length(Q_succ)
                    VN(i)=Wp(Q_succ(i)); 
                end
                Vmax(p)=max(VN);    
             else
                Vmax(p)=2*Hmax;  
            end            
        end
        [v,c]=min(Vmax);
        Wn(q)=min(Vn(q),v); 
        
        if Wn(q)~=Wp(q)
            K(q)=iter;
            C_mF(q,:)=(Vmax==v);
        end
    end             
    fixpoint=isequal(Wp,Wn);
end
disp('Fixed point reached')
Iterations=iter


toc
           
disp('Abstraction and controller synthesis without Nested approach finished')
%load handel
%sound(y,Fs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots: attractor and basin of attraction without partition refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Q0=(Vn==min(Vn)); %Attractor
Attract=find(Q0);
Q2=(Wn==min(Vn)); %Basin of attraction
Q=(Vn==min(Vn))+(Wn==min(Vn));

QVp=find(Vn==min(Vn));
QWp=find(Wn==min(Vn));

Q_flat=zeros(n_x(1),n_x(2),n_x(3));
for q=1:n_states
    c_q=state2coord(q,p_x);
    Q_flat(c_q(1),c_q(2),c_q(3))=Q(q);
end

%%%%%%%
% 3 Dimension
%%%%%%%
figure(12) % Basin of attraction
for i=1:max(size(QWp))
    c_q=state2coord(QWp(i),p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    cube_plot((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'yellow')
end
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2) bound(3,1) bound(3,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
zlabel('x_3','FontSize',18);
axis([-4 4 -3 3 -4 4]);


figure(11) % Attractor
for i=1:max(size(QVp))
    c_q=state2coord(QVp(i),p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    cube_plot((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'green')
end
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2) bound(3,1) bound(3,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
zlabel('x_3','FontSize',18);

%%%%%%%
% 2 Dimension
%%%%%%%

t1=floor((0-bound(3,1))/d_x(3))+1;

figure(13)
title('x_3=0','Fontsize',18)
for i=1:n_states
    c_q=state2coord(i,p_x);
    if c_q(3,1)==t1
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'white')
    end
end
hold on
for i=1:max(size(QWp))
    c_q=state2coord(QWp(i),p_x);
    if c_q(3,1)==t1
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'yellow')
    end
end
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
hold on
for i=1:max(size(Attract))
    c_q=state2coord(Attract(i),p_x);
    if c_q(3,1)==t1
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'green')
    end
end
%obstacles
Rectangle_plot([obstacleState_1(1,1);obstacleState_1(2,1)],obstacleState_1(1,2)-obstacleState_1(1,1),obstacleState_1(2,2)-obstacleState_1(2,1),[128, 64, 0]*1/255)
Rectangle_plot([obstacleState_2(1,1);obstacleState_2(2,1)],obstacleState_2(1,2)-obstacleState_2(1,1),obstacleState_2(2,2)-obstacleState_2(2,1),[128, 64, 0]*1/255)
%border of the target set
plot([X_refSet(1,1) X_refSet(1,2) X_refSet(1,2) X_refSet(1,1) X_refSet(1,1)],[X_refSet(2,1) X_refSet(2,1) X_refSet(2,2) X_refSet(2,2) X_refSet(2,1)],'b','LineWidth',2)
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the abstraction for the nested transition relation which depends of
% the previous controller matrix C_mF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial space discretization
r=4;                              % parameter for the new abstraction
n_xN=r*n_x;     % number of cells in each direction
p_xN=[1;cumprod(n_xN)];             % [1;n_x(1);n_x(1)*n_x(2)...] used for coordinate/state correspondance 
n_statesN=p_xN(end);                % total number of states
d_xN=(bound(:,2)-bound(:,1))./n_xN; 


Q0=(Vn==min(Vn));
Attract=find(Q0);

tic

disp('Step 1: Computation of the abstraction''s nested transition relation')
h = waitbar(0,'Abstraction''s nested in progress...');

Delta_N=zeros(n_statesN*n_modes,2);  % Delta contains the transition relation
H_N=zeros(n_statesN,1); % H contains distance to safe set

for q=1:n_statesN 
    c_q=state2coord(q,p_xN); % coordinates of the cell
    x_c=bound(:,1)+(c_q-0.5).*d_xN; % center of rectangle associated to state q %nx1 (in the actual state space)
    H_N(q)=distanceVessel(x_c,0.5*d_xN,bound,X_ref,obstacleState_1,obstacleState_2);   
end

NC_N=[];

for q=1:max(size(Attract))
    waitbar(q /max(size(Attract)),h)
    c_q=state2coord(Attract(q),p_x);
    modes_N=find(C_mF(Attract(q),:));
    a=r*(c_q(1)-1)+1;
    b=r*(c_q(2)-1)+1;
    c=r*(c_q(3)-1)+1;
    [A1, B1, C1]=ndgrid(a:a+r-1,b:b+r-1,c:c+r-1);
    d=[A1(:),B1(:),C1(:)]';
    for i=1:max(size(d(1,:)))
        c_qN=d(:,i);
        qN=coord2state(c_qN,p_xN);
        x_c=bound(:,1)+(c_qN-0.5).*d_xN; 
        for p=1:length(modes_N)
            u_l=mode(:,modes_N(p));
            x_succ_c=F(x_c, u_l);
            d_x_succ=0.5*A(u_l)*d_xN;
            x_succ=[x_succ_c-d_x_succ,x_succ_c+d_x_succ]; 
            if x_succ(3,1)+tol<-pi 
                x_succ(3,1)=x_succ(3,1)+2*pi;
            end
            if x_succ(3,2)+tol<-pi       
                x_succ(3,2)=x_succ(3,2)+2*pi;
            end
            if x_succ(3,2)-tol>pi 
                x_succ(3,2)=x_succ(3,2)-2*pi;
            end
            if x_succ(3,1)-tol>pi        
                x_succ(3,1)=x_succ(3,1)-2*pi;
            end
            
            c_ql=state2coord(Delta((Attract(q)-1)*n_modes+modes_N(p),1),p_x);
            c_qu=state2coord(Delta((Attract(q)-1)*n_modes+modes_N(p),2),p_x);
            x_l=(bound(:,1)+(c_ql-1).*d_x); 
            x_u=(bound(:,1)+(c_qu).*d_x);
            if max(x_l>(x_succ(:,1)+tol)) || max((x_succ(:,2)-tol)>x_u)
                NC_N=[NC_N;q];
            end    
            
            if bound(1,1)<=(x_succ(1,1)+tol) && (x_succ(1,2)-tol)<=bound(1,2) && bound(2,1)<=(x_succ(2,1)+tol) && (x_succ(2,2)-tol)<=bound(2,2)...
                && ((x_succ(1,2)-tol)<=obstacleState_1(1,1) || (x_succ(2,2)-tol)<=obstacleState_1(2,1) || (x_succ(1,1)+tol)>=obstacleState_1(1,2) || (x_succ(2,1)+tol)>=obstacleState_1(2,2))...
                && ((x_c(1)+0.5*d_x(1)-tol)<=obstacleState_1(1,1) || (x_c(2)+0.5*d_x(2)-tol)<=obstacleState_1(2,1) || (x_c(1)-0.5*d_x(1)+tol)>=obstacleState_1(1,2) || (x_c(2)-0.5*d_x(1)+tol)>=obstacleState_1(2,2))...
                && ((x_succ(1,2)-tol)<=obstacleState_2(1,1) || (x_succ(2,2)-tol)<=obstacleState_2(2,1) || (x_succ(1,1)+tol)>=obstacleState_2(1,2) || (x_succ(2,1)+tol)>=obstacleState_2(2,2))...
                && ((x_c(1)+0.5*d_x(1)-tol)<=obstacleState_2(1,1) || (x_c(2)+0.5*d_x(2)-tol)<=obstacleState_2(2,1) || (x_c(1)-0.5*d_x(1)+tol)>=obstacleState_2(1,2) || (x_c(2)-0.5*d_x(1)+tol)>=obstacleState_2(2,2))
                min_succ=coord2state(floor((x_succ(:,1)-bound(:,1))./d_xN)+1,p_xN); 
                max_succ=coord2state(ceil((x_succ(:,2)-bound(:,1))./d_xN),p_xN);
                Delta_N((qN-1)*n_modes+modes_N(p),:)=[min_succ,max_succ];
            end
        end   
    end     
end
close(h) 
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantitative synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic


% Safety synthesis

Hmax_N=max(H_N); % Maximal distance to safe set with the partition refinement
Vn_N=H_N;
Vp_N=zeros(n_statesN,1);
C_mF_N=zeros(n_statesN,n_modes); % controller matrix for the new controller
iter=0;

disp('Step 2: Safety Controller synthesis for nested system')
disp('Fixed point algorithm')
fixpoint=0;
while ~fixpoint
    iter=iter+1
    Vp_N=Vn_N;
    for q=1:max(size(Attract))
        c_q=state2coord(Attract(q),p_x);
        a=r*(c_q(1)-1)+1;
        b=r*(c_q(2)-1)+1;
        c=r*(c_q(3)-1)+1;
        [A1, B1, C1]=ndgrid(a:a+r-1,b:b+r-1,c:c+r-1);
        d=[A1(:),B1(:),C1(:)]';
        modes_N=find(C_mF(Attract(q),:));
        
        for i=1:max(size(d(1,:)))
            c_qN=d(:,i);
            qN=coord2state(c_qN,p_xN);    
            Vmax=ones(1,n_modes)*2*Hmax_N; 
            for p=1:length(modes_N)
                if any(Delta_N((qN-1)*n_modes+modes_N(p),:)) 
                    Imin=state2coord(Delta_N((qN-1)*n_modes+modes_N(p),1),p_xN);
                    Imax=state2coord(Delta_N((qN-1)*n_modes+modes_N(p),2),p_xN);
                    if Imin(3)<=Imax(3)
                        Q_succ=coordbox2state(Imin,Imax,p_xN)';             
                    else
                        Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_xN)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_xN)'];
                    end
                
                    VN=[]; 
                    for i=1:length(Q_succ)
                         VN(i)=Vp_N(Q_succ(i)); 
                    end
                        Vmax(modes_N(p))=max(VN); 
                end            
             end
             [v,c]=min(Vmax);
             Vn_N(qN)=max(H_N(qN),v); 
         end
    end           
    fixpoint=isequal(Vp_N,Vn_N);
end

%Finding all available controllers for safety specification
for q=1:max(size(Attract))
        modes_N=find(C_mF(Attract(q),:));
        c_q=state2coord(Attract(q),p_x);
        a=r*(c_q(1)-1)+1;
        b=r*(c_q(2)-1)+1;
        c=r*(c_q(3)-1)+1;
        [A1, B1,C1]=ndgrid(a:a+r-1,b:b+r-1,c:c+r-1);
        d=[A1(:),B1(:),C1(:)]';
        
        for i=1:max(size(d(1,:)))
            c_qN=d(:,i);
            qN=coord2state(c_qN,p_xN);
            for p=1:length(modes_N)
                 if any(Delta_N((qN-1)*n_modes+modes_N(p),:)) 
                    Imin=state2coord(Delta_N((qN-1)*n_modes+modes_N(p),1),p_xN);
                    Imax=state2coord(Delta_N((qN-1)*n_modes+modes_N(p),2),p_xN);
                    if Imin(3)<=Imax(3)
                        Q_succ=coordbox2state(Imin,Imax,p_xN)';             
                    else
                        Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_xN)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_xN)'];
                    end
                      VN=[]; 
                      for i=1:length(Q_succ)
                           VN(i)=Vn_N(Q_succ(i)); 
                      end
                      if max(VN)<=Vn_N(qN) 
                            C_mF_N(qN,modes_N(p))=1; 
                      end    
                 end            
            end
        end
end    

disp('Fixed point reached')
Iterations=iter

% Reachability synthesis

% Value function
Wn_N=Vn_N;
Wp_N=zeros(n_statesN,1);
K_N=zeros(n_statesN,1); 

iter=0;
disp('Step 3: Reachability Controller synthesis for nested system')
disp('Fixed point algorithm')
fixpoint=0;
while ~fixpoint
     iter=iter+1
     Wp_N=Wn_N;       
    for q=1:max(size(Attract))
        c_q=state2coord(Attract(q),p_x);
        a=r*(c_q(1)-1)+1;
        b=r*(c_q(2)-1)+1;
        c=r*(c_q(3)-1)+1;
        [A1, B1, C1]=ndgrid(a:a+r-1,b:b+r-1,c:c+r-1);
        d=[A1(:),B1(:),C1(:)]';
        modes_N=find(C_mF(Attract(q),:));
        
        for i=1:max(size(d(1,:)))
            c_qN=d(:,i);
            qN=coord2state(c_qN,p_xN);
            Vmax=ones(1,n_modes)*2*Hmax_N; 
            for p=1:length(modes_N)
            	if any(Delta_N((qN-1)*n_modes+modes_N(p),:)) 
                    Imin=state2coord(Delta_N((qN-1)*n_modes+modes_N(p),1),p_xN);
                    Imax=state2coord(Delta_N((qN-1)*n_modes+modes_N(p),2),p_xN);
                    if Imin(3)<=Imax(3)
                        Q_succ=coordbox2state(Imin,Imax,p_xN)';             
                    else
                        Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_xN)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_xN)'];
                    end   
                    VN=[]; 
                    for i=1:length(Q_succ)
                          VN(i)=Wp_N(Q_succ(i)); 
                    end
                    Vmax(modes_N(p))=max(VN);   
                else
                    Vmax(modes_N(p))=2*Hmax_N;  
                end            
             end
             [v,c]=min(Vmax);
             Wn_N(qN)=min(Vn_N(qN),v); 
             if Wn_N(qN)~=Wp_N(qN)
             	K_N(qN)=iter;           
                C_mF_N(qN,:)=(Vmax==v);
             end
        end
    end   
    fixpoint=isequal(Wp_N,Wn_N);
end

disp('Fixed point reached')
Iterations=iter
% 
toc


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots: basin of attraction and attractors with partition refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%﻿First Attractor in new partition
VnNR=[];
WnNR=[];
for q=1:max(size(Attract))
    c_q=state2coord(Attract(q),p_x);
    a=r*(c_q(1)-1)+1;
    b=r*(c_q(2)-1)+1;
    c=r*(c_q(3)-1)+1;
    [A1, B1,C1]=ndgrid(a:a+r-1,b:b+r-1,c:c+r-1);
    d=[A1(:),B1(:),C1(:)]';
    for i=1:max(size(d(1,:)))
        c_qN=d(:,i);
        qN=coord2state(c_qN,p_xN);
        VnNR=[VnNR;qN Vn_N(qN)];
        WnNR=[WnNR;qN Wn_N(qN)];
    end
end

%%% New Attractor in new a partition
QVp_N1=find(VnNR(:,2)-tol<=max(WnNR(:,2)));
QVp_N2=find(WnNR(:,2)-tol<=max(WnNR(:,2)));


%%% Plot 3D

figure(20) % Basin of attraction
for i=1:n_states
    c_q=state2coord(i,p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    cube_plotTrans2((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'white')
end
for i=1:max(size(QWp))
    c_q=state2coord(QWp(i),p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    cube_plot((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'yellow')
end
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2) bound(3,1) bound(3,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
zlabel('x_3','FontSize',18);

figure(21) %first attractor in the partition refinement
for i=1:n_states
    c_q=state2coord(i,p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    cube_plotTrans2((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'white')
end
for i=1:max(size(QVp_N2))
    c_q=state2coord(WnNR(QVp_N2(i),1),p_xN);
    x_l=bound(:,1)+(c_q-0.5).*d_xN;
    cube_plot((x_l-0.5*d_xN)',d_xN(1),d_xN(2),d_xN(3),'green')
end
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2) bound(3,1) bound(3,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
zlabel('x_3','FontSize',18);

figure(22)
for i=1:n_states
    c_q=state2coord(i,p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    cube_plotTrans2((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'white')
end
for i=1:max(size(QVp_N1))
    c_q=state2coord(VnNR(QVp_N1(i),1),p_xN);
    x_l=bound(:,1)+(c_q-0.5).*d_xN;
    cube_plot((x_l-0.5*d_xN)',d_xN(1),d_xN(2),d_xN(3),'red')
end
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2) bound(3,1) bound(3,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
zlabel('x_3','FontSize',18);

figure(23)
for i=1:n_states
    c_q=state2coord(i,p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    cube_plotTrans2((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'white')
end
for i=1:max(size(QWp))
    c_q=state2coord(QWp(i),p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    cube_plotTrans1((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'yellow')
end
for i=1:max(size(QVp_N2))
    c_q=state2coord(WnNR(QVp_N2(i),1),p_xN);
    x_l=bound(:,1)+(c_q-0.5).*d_xN;
    cube_plotTrans1((x_l-0.5*d_xN)',d_xN(1),d_xN(2),d_xN(3),'green')
end
for i=1:max(size(QVp_N1))
    c_q=state2coord(VnNR(QVp_N1(i),1),p_xN);
    x_l=bound(:,1)+(c_q-0.5).*d_xN;
    cube_plot((x_l-0.5*d_xN)',d_xN(1),d_xN(2),d_xN(3),'red')
end
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2) bound(3,1) bound(3,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
zlabel('x_3','FontSize',18);


figure(230)
mid=(X_refSet(1:2,1)+X_refSet(1:2,2))/2;
for i=1:n_states
    c_q=state2coord(i,p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    if mid(1,1)<x_l(1,1) || mid(2,1)<x_l(2,1)
    cube_plotTrans2((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'white')
    end
end
for i=1:max(size(QWp))
    c_q=state2coord(QWp(i),p_x);
    x_l=bound(:,1)+(c_q-0.5).*d_x;
    if mid(1,1)<x_l(1,1) || mid(2,1)<x_l(2,1)
    cube_plotTrans1((x_l-0.5*d_x)',d_x(1),d_x(2),d_x(3),'yellow')
    end
end
for i=1:max(size(QVp_N2))
    c_q=state2coord(WnNR(QVp_N2(i),1),p_xN);
    x_l=bound(:,1)+(c_q-0.5).*d_xN;
    if mid(1,1)<x_l(1,1) || mid(2,1)<x_l(2,1)
    cube_plotTrans1((x_l-0.5*d_xN)',d_xN(1),d_xN(2),d_xN(3),'green')
    end
end
for i=1:max(size(QVp_N1))
    c_q=state2coord(VnNR(QVp_N1(i),1),p_xN);
    x_l=bound(:,1)+(c_q-0.5).*d_xN;
    if mid(1,1)<x_l(1,1) || mid(2,1)<x_l(2,1)
    cube_plot((x_l-0.5*d_xN)',d_xN(1),d_xN(2),d_xN(3),'red')
    end
end
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2) bound(3,1) bound(3,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
zlabel('x_3','FontSize',18);

%%%%%%%
% 2 Dimension: cuts at particular values of x3
%%%%%%%

t1=floor(((-pi/7)-bound(3,1))/d_x(3))+1;
t01=floor(((-pi/7)-bound(3,1))/d_xN(3))+1;

t2=floor((0-bound(3,1))/d_x(3))+1;
t02=floor((0-bound(3,1))/d_xN(3))+1;

t3=floor(((pi/10)-bound(3,1))/d_x(3))+1;
t03=floor(((pi/10)-bound(3,1))/d_xN(3))+1;


figure(24)
title('x_3=-\pi/7','Fontsize',18)
for i=1:n_states
    c_q=state2coord(i,p_x);
    if c_q(3,1)==t1
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'white')
    end
end
hold on
for i=1:max(size(QWp))
    c_q=state2coord(QWp(i),p_x);
    if c_q(3,1)==t1
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'yellow')
    end
end
for p=1:max(size(Attract))
    c_q=state2coord(Attract(p),p_x);
    a=r*(c_q(1)-1)+1;
    b=r*(c_q(2)-1)+1;
    c=r*(c_q(3)-1)+1;
    [A1, B1, C1]=ndgrid(a:a+r-1,b:b+r-1,c:c+r-1);
    d=[A1(:),B1(:),C1(:)]';  
    for i=1:max(size(d(1,:)))
        c_qN=d(:,i);
        qN=coord2state(c_qN,p_xN);
        if c_qN(3,1)==t01
            x_l=bound(:,1)+(c_qN-0.5).*d_xN;
            xx=x_l(1:2,1);
            Rectangle_plot((xx-0.5*d_xN(1:2,1))',d_xN(1),d_xN(2),'green')
        end
    end
end  
for i=1:max(size(QVp_N1))
    c_qN=state2coord(VnNR(QVp_N1(i),1),p_xN);
        if c_qN(3,1)==t01
            x_l=bound(:,1)+(c_qN-0.5).*d_xN;
            xx=x_l(1:2,1);
            Rectangle_plot((xx-0.5*d_xN(1:2,1))',d_xN(1),d_xN(2),'red')
        end
end
Rectangle_plot([obstacleState_1(1,1);obstacleState_1(2,1)],obstacleState_1(1,2)-obstacleState_1(1,1),obstacleState_1(2,2)-obstacleState_1(2,1),[128, 64, 0]*1/255)
Rectangle_plot([obstacleState_2(1,1);obstacleState_2(2,1)],obstacleState_2(1,2)-obstacleState_2(1,1),obstacleState_2(2,2)-obstacleState_2(2,1),[128, 64, 0]*1/255)
plot([X_refSet(1,1) X_refSet(1,2) X_refSet(1,2) X_refSet(1,1) X_refSet(1,1)],[X_refSet(2,1) X_refSet(2,1) X_refSet(2,2) X_refSet(2,2) X_refSet(2,1)],'b','LineWidth',2)
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);

figure(25)
title('x_3=0','Fontsize',18)
for i=1:n_states
    c_q=state2coord(i,p_x);
    if c_q(3,1)==t1
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'white')
    end
end
hold on
for i=1:max(size(QWp))
    c_q=state2coord(QWp(i),p_x);
    if c_q(3,1)==t2
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'yellow')
    end
end  
for p=1:max(size(Attract))
    c_q=state2coord(Attract(p),p_x);
    a=r*(c_q(1)-1)+1;
    b=r*(c_q(2)-1)+1;
    c=r*(c_q(3)-1)+1;
    [A1, B1, C1]=ndgrid(a:a+r-1,b:b+r-1,c:c+r-1);
    d=[A1(:),B1(:),C1(:)]';  
    for i=1:max(size(d(1,:)))
        c_qN=d(:,i);
        qN=coord2state(c_qN,p_xN);
        if c_qN(3,1)==t02
            x_l=bound(:,1)+(c_qN-0.5).*d_xN;
            xx=x_l(1:2,1);
            Rectangle_plot((xx-0.5*d_xN(1:2,1))',d_xN(1),d_xN(2),'green')
        end
    end
end
for i=1:max(size(QVp_N1))
    c_qN=state2coord(VnNR(QVp_N1(i),1),p_xN);
        if c_qN(3,1)==t02
            x_l=bound(:,1)+(c_qN-0.5).*d_xN;
            xx=x_l(1:2,1);
            Rectangle_plot((xx-0.5*d_xN(1:2,1))',d_xN(1),d_xN(2),'red')
        end
end
Rectangle_plot([obstacleState_1(1,1);obstacleState_1(2,1)],obstacleState_1(1,2)-obstacleState_1(1,1),obstacleState_1(2,2)-obstacleState_1(2,1),[128, 64, 0]*1/255)
Rectangle_plot([obstacleState_2(1,1);obstacleState_2(2,1)],obstacleState_2(1,2)-obstacleState_2(1,1),obstacleState_2(2,2)-obstacleState_2(2,1),[128, 64, 0]*1/255)
plot([X_refSet(1,1) X_refSet(1,2) X_refSet(1,2) X_refSet(1,1) X_refSet(1,1)],[X_refSet(2,1) X_refSet(2,1) X_refSet(2,2) X_refSet(2,2) X_refSet(2,1)],'b','LineWidth',2)
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);

figure(26)
title('x_3=\pi/10','Fontsize',18)
for i=1:n_states
    c_q=state2coord(i,p_x);
    if c_q(3,1)==t1
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'white')
    end
end
hold on
for i=1:max(size(QWp))
    c_q=state2coord(QWp(i),p_x);
    if c_q(3,1)==t3
       x_l=bound(:,1)+(c_q-0.5).*d_x;
       xx=x_l(1:2,1);
       Rectangle_plot((xx-0.5*d_x(1:2,1))',d_x(1),d_x(2),'yellow')
    end
end
for p=1:max(size(Attract))
    c_q=state2coord(Attract(p),p_x);
    a=r*(c_q(1)-1)+1;
    b=r*(c_q(2)-1)+1;
    c=r*(c_q(3)-1)+1;
    [A1, B1, C1]=ndgrid(a:a+r-1,b:b+r-1,c:c+r-1);
    d=[A1(:),B1(:),C1(:)]';  
    for i=1:max(size(d(1,:)))
        c_qN=d(:,i);
        qN=coord2state(c_qN,p_xN);
        if c_qN(3,1)==t03
            x_l=bound(:,1)+(c_qN-0.5).*d_xN;
            xx=x_l(1:2,1);
            Rectangle_plot((xx-0.5*d_xN(1:2,1))',d_xN(1),d_xN(2),'green')
        end
    end
end
for i=1:max(size(QVp_N1))
    c_qN=state2coord(VnNR(QVp_N1(i),1),p_xN);
        if c_qN(3,1)==t03
            x_l=bound(:,1)+(c_qN-0.5).*d_xN;
            xx=x_l(1:2,1);
            Rectangle_plot((xx-0.5*d_xN(1:2,1))',d_xN(1),d_xN(2),'red')
        end
end
Rectangle_plot([obstacleState_1(1,1);obstacleState_1(2,1)],obstacleState_1(1,2)-obstacleState_1(1,1),obstacleState_1(2,2)-obstacleState_1(2,1),[128, 64, 0]*1/255)
Rectangle_plot([obstacleState_2(1,1);obstacleState_2(2,1)],obstacleState_2(1,2)-obstacleState_2(1,1),obstacleState_2(2,2)-obstacleState_2(2,1),[128, 64, 0]*1/255)
plot([X_refSet(1,1) X_refSet(1,2) X_refSet(1,2) X_refSet(1,1) X_refSet(1,1)],[X_refSet(2,1) X_refSet(2,1) X_refSet(2,2) X_refSet(2,2) X_refSet(2,1)],'b','LineWidth',2)
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);


%%%% Trajectories in the state space

x=[2;2;pi/3];
traj=[x];
N=60;
TV=[];
Tmode=[];
for t=1:N 
   q=coord2state(floor((x-bound(:,1))./d_x)+1,p_x);
   if ismember(q,QVp)
       qN=coord2state(floor((x-bound(:,1))./d_xN)+1,p_xN);
       m=find(C_mF_N(qN,:),1,'last');
   elseif ismember(q,QWp)
       m=find(C_mF(q,:),1,'last');
   else
       disp('Reached uncontrollable state')
       break
   end
  
   x=F(x,mode(:,m));
   if x(3)>pi
       x(3)=x(3)-2*pi;
   elseif x(3)<-pi
       x(3)=x(3)+2*pi;
   end
   Tmode=[Tmode,m];
   traj=[traj x];
end

figure(27)
plot(traj(1,:),traj(2,:),'black','LineWidth',2)
hold on
Rectangle_plot([obstacleState_1(1,1);obstacleState_1(2,1)],obstacleState_1(1,2)-obstacleState_1(1,1),obstacleState_1(2,2)-obstacleState_1(2,1),[128, 64, 0]*1/255)
Rectangle_plot([obstacleState_2(1,1);obstacleState_2(2,1)],obstacleState_2(1,2)-obstacleState_2(1,1),obstacleState_2(2,2)-obstacleState_2(2,1),[128, 64, 0]*1/255)
plot([X_refSet(1,1) X_refSet(1,2) X_refSet(1,2) X_refSet(1,1) X_refSet(1,1)],[X_refSet(2,1) X_refSet(2,1) X_refSet(2,2) X_refSet(2,2) X_refSet(2,1)],'b','LineWidth',2)
plot(traj(1,1),traj(2,1),'o','color','black','LineWidth',15)
axis([bound(1,1) bound(1,2) bound(2,1) bound(2,2)])
set(gca,'FontSize',18)
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);

%%%% Trajectories with respect to the time
figure3 = figure(28);
% Create subplot
subplot1 = subplot(5,1,1,'Parent',figure3,'FontSize',14);
ylim([-pi pi])
box(subplot1,'on');
hold(subplot1,'on');
% Create plot
plot([0:t]*T0,traj(3,:),'Parent',subplot1,'Color',[0 0 1],'LineWidth',2);
plot([0:t]*T0,0*traj(1,:) + X_refSet(3,1),'-.','Parent',subplot1,'Color','red','LineWidth',1);
plot([0:t]*T0,0*traj(1,:) + X_refSet(3,2),'-.','Parent',subplot1,'Color','red','LineWidth',1);
ylabel('x_3','FontSize',14);
set(gca,'XTick',[])
set(gca,'FontSize',18)

% Create subplot
subplot1 = subplot(5,1,2,'Parent',figure3,'FontSize',14);
ylim([bound(1,1) bound(1,2)])
box(subplot1,'on');
hold(subplot1,'on');
plot([0:t]*T0,traj(1,:),'Parent',subplot1,'Color',[0 0 1],'LineWidth',2);
plot([0:t]*T0,traj(2,:),'Parent',subplot1,'Color',[127,179,131]*1/250,'LineWidth',2);
plot([0:t]*T0,0*traj(1,:) + X_refSet(1,1),'-.','Parent',subplot1,'Color','red','LineWidth',1);
plot([0:t]*T0,0*traj(1,:) + X_refSet(1,2),'-.','Parent',subplot1,'Color','red','LineWidth',1);
ylabel('x_1, x_2','FontSize',14);
set(gca,'XTick',[])
set(gca,'FontSize',18)
% Create subplot
subplot1 = subplot(5,1,3,'Parent',figure3,'FontSize',14);
ylim([inputbound(1,1) inputbound(1,2)])
box(subplot1,'on');
hold(subplot1,'on');
plot([0:t-1]*T0,mode(1,Tmode),'Parent',subplot1,'Color',[0 0 1],'LineWidth',2);
ylabel('u_1','FontSize',14);
set(gca,'XTick',[])
set(gca,'FontSize',18)

% Create subplot
subplot1 = subplot(5,1,4,'Parent',figure3,'FontSize',14);
ylim([inputbound(2,1) inputbound(2,2)])
box(subplot1,'on');
hold(subplot1,'on');
plot([0:t-1]*T0,mode(2,Tmode),'Parent',subplot1,'Color',[0 0 1],'LineWidth',2);
ylabel('u_2','FontSize',14);
set(gca,'XTick',[])
set(gca,'FontSize',18)
 
% Create subplot
subplot1 = subplot(5,1,5,'Parent',figure3,'FontSize',14);
ylim([inputbound(3,1) inputbound(3,2)])
box(subplot1,'on');
hold(subplot1,'on');
plot([0:t-1]*T0,mode(3,Tmode),'Parent',subplot1,'Color',[0 0 1],'LineWidth',2);
xlabel('t (min)','FontSize',14);
ylabel('u_3','FontSize',14);
set(gca,'FontSize',18)

%%
%%%% This is not necessary, this is only to get a performance comparison
%%%% with respect to the computation time
% Abstraction and control synthesis with the fine uniform partition  

% Quantitative Reachability synthesis
clear
close all
clc
tol=10^-6; %tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Bounds
bound=[-3.5 3.5;-2.6 2.6;-pi pi];
inputbound=[-0.7 2;-1 1;-1 1];
disturbound=[-0.001 0.001;-0.001 0.001;-0.001 0.001];


obstacleState_1=[0.7 1;0.7 2.6;-pi pi];
obstacleState_2=[-1.8 -1.5;-2.6 -1;-pi pi];


% Target set

X_ref=[-0.5;-0.5;0];
X_refSet=[-1.2 0.1;-1.2 0.1;-0.5 0.5];
% Dynamics
%tol=10^-6;
T0 =1;

F = @(x, u) [x(1)+T0*(u(1)*cos(x(3))-u(2)*sin(x(3)));...
             x(2)+T0*(u(1)*sin(x(3))+u(2)*cos(x(3)));...
             x(3)+T0*u(3)];
         
A = @(u) [1 0 T0*(abs(u(1))+abs(u(2)));...
          0 1 T0*(abs(u(1))+abs(u(2)));...
          0 0 1];
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Abstraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% Space discretization
r=4;
n_x=[15;11;10]*r; % to get the fine partition refinement
p_x=[1;cumprod(n_x)];
n_states=p_x(end);                % total number of states
d_x=(bound(:,2)-bound(:,1))./n_x;    

% Input discretization
n_modes=[6;5;5];       % number of modes
d_u=(inputbound(:,2)-inputbound(:,1))./(n_modes-1);
mode=[];
for i=1:n_modes(1)
    for j=1:n_modes(2)
        for k=1:n_modes(3)
            mode=[mode inputbound(:,1)+[i-1;j-1;k-1].*d_u];
        end    
    end
end
n_modes=length(mode);

disp('Step 1: Computation of the abstraction''s transition relation')
h = waitbar(0,'Abstraction in progress...');

H=zeros(n_states,1); % H contains distance to safe set
Delta=zeros(n_states*n_modes,2);  % Delta contains the transition relation

for q=1:n_states 
    waitbar(q / n_states,h)
    c_q=state2coord(q,p_x); % coordinates of the cell
    x_c=bound(:,1)+(c_q-0.5).*d_x; % center of rectangle associated to state q %2x1 (in the actual state space)
    % Compute the distance to the safe set
    H(q)=distanceVessel(x_c,0.5*d_x,bound,X_ref,obstacleState_1,obstacleState_2);
    % Compute the transition relation
    for p=1:n_modes
        x_succ_c=F(x_c, mode(:,p));
        d_x_succ=0.5*A(mode(:,p))*d_x;
        x_succ=[x_succ_c-d_x_succ,x_succ_c+d_x_succ]; % Over-approximation of the reachable set from x under mode p
        if x_succ(3,1)-tol<-pi 
            x_succ(3,1)=x_succ(3,1)+2*pi;
        end
        if x_succ(3,2)-tol<-pi       %to cover all posibilities 
            x_succ(3,2)=x_succ(3,2)+2*pi;
        end
        if x_succ(3,2)+tol>pi 
            x_succ(3,2)=x_succ(3,2)-2*pi;
        end
        if x_succ(3,1)+tol>pi        %to cover all posibilities
            x_succ(3,1)=x_succ(3,1)-2*pi;
        end
        
         if bound(1,1)<=(x_succ(1,1)+tol) && (x_succ(1,2)-tol)<=bound(1,2) && bound(2,1)<=(x_succ(2,1)+tol) && (x_succ(2,2)-tol)<=bound(2,2)...
            && ((x_succ(1,2)-tol)<=obstacleState_1(1,1) || (x_succ(2,2)-tol)<=obstacleState_1(2,1) || (x_succ(1,1)+tol)>=obstacleState_1(1,2) || (x_succ(2,1)+tol)>=obstacleState_1(2,2))...
            && ((x_c(1)+0.5*d_x(1)-tol)<=obstacleState_1(1,1) || (x_c(2)+0.5*d_x(2)-tol)<=obstacleState_1(2,1) || (x_c(1)-0.5*d_x(1)+tol)>=obstacleState_1(1,2) || (x_c(2)-0.5*d_x(1)+tol)>=obstacleState_1(2,2))...
            && ((x_succ(1,2)-tol)<=obstacleState_2(1,1) || (x_succ(2,2)-tol)<=obstacleState_2(2,1) || (x_succ(1,1)+tol)>=obstacleState_2(1,2) || (x_succ(2,1)+tol)>=obstacleState_2(2,2))...
            && ((x_c(1)+0.5*d_x(1)-tol)<=obstacleState_2(1,1) || (x_c(2)+0.5*d_x(2)-tol)<=obstacleState_2(2,1) || (x_c(1)-0.5*d_x(1)+tol)>=obstacleState_2(1,2) || (x_c(2)-0.5*d_x(1)+tol)>=obstacleState_2(2,2))
            min_succ=coord2state(floor((x_succ(:,1)-bound(:,1))./d_x)+1,p_x); 
            max_succ=coord2state(ceil((x_succ(:,2)-bound(:,1))./d_x),p_x);
            Delta((q-1)*n_modes+p,:)=[min_succ,max_succ];   
        end
    end
end
close(h) 
toc

NC=[];
for q=1:n_states
         if max(~any(Delta((q-1)*n_modes+1:(q-1)*n_modes+n_modes,1)))
             NC=[NC;q];
         end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantitative synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
Hmax=max(H);%Inf; % Maximal distance to safe set

% Safety synthesis
Vn=H;
Vp=zeros(n_states,1);
C_mF=zeros(n_states,n_modes);
C=zeros(n_states,1);
Vmax=zeros(1,n_modes);
iter=0;

disp('Step 2: Safety Controller synthesis')
disp('Fixed point algorithm')
fixpoint=0;
while ~fixpoint
    iter=iter+1
    Vp=Vn;
    for q=1:n_states
        Vmax=zeros(1,n_modes);
        for p=1:n_modes
            if any(Delta((q-1)*n_modes+p,:)) % then we have successors
                Imin=state2coord(Delta((q-1)*n_modes+p,1),p_x);
                Imax=state2coord(Delta((q-1)*n_modes+p,2),p_x);
                if Imin(3)<=Imax(3)
                    Q_succ=coordbox2state(Imin,Imax,p_x)';             
                else
                    Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_x)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_x)'];
                end
                VN=[]; % This was the bug! need to clear this variable, because of the different lengths!
                for i=1:length(Q_succ)
                    VN(i)=Vp(Q_succ(i));
                end
                Vmax(p)=max(VN);  % Max value at successor states  
             else
                Vmax(p)=2*Hmax;
            end
        end
        [v,c]=min(Vmax);
        Vn(q)=max(H(q),v);  % MinMax value over inputs %NEW
        if any(Delta((q-1)*n_modes+c,:))
          C(q)=c;
        else
          C(q)=0;
        end
    end             
    fixpoint=isequal(Vp,Vn);
end

%Finding the controller for safety
for q=1:n_states
        for p=1:n_modes
            if any(Delta((q-1)*n_modes+p,:)) 
                Imin=state2coord(Delta((q-1)*n_modes+p,1),p_x);
                Imax=state2coord(Delta((q-1)*n_modes+p,2),p_x);
                if Imin(3)<=Imax(3)
                    Q_succ=coordbox2state(Imin,Imax,p_x)';             
                else
                    Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_x)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_x)'];
                end  
                VN=[]; 
                for i=1:length(Q_succ)
                    VN(i)=Vn(Q_succ(i)); 
                end
                if max(VN)<=Vn(q)  
                   C_mF(q,p)=1;  
                end    
            end     
        end
end


disp('Fixed point reached')
Iterations=iter

% Reachability synthesis
Wn=Vn;
Wp=zeros(n_states,1);
K=zeros(n_states,1);

iter=0;

disp('Step 3: Reachability Controller synthesis')
disp('Fixed point algorithm')
fixpoint=0;
while ~fixpoint
    iter=iter+1
    Wp=Wn;
    for q=1:n_states
        for p=1:n_modes
            if any(Delta((q-1)*n_modes+p,:)) 
                Imin=state2coord(Delta((q-1)*n_modes+p,1),p_x);
                Imax=state2coord(Delta((q-1)*n_modes+p,2),p_x);
                if Imin(3)<=Imax(3)
                    Q_succ=coordbox2state(Imin,Imax,p_x)';             
                else
                    Q_succ=[coordbox2state(Imin,[Imax(1);Imax(2);n_x(3)],p_x)',...
                            coordbox2state([Imin(1);Imin(2);1],Imax,p_x)'];
                end
                VN=[]; 
                for i=1:length(Q_succ)
                    VN(i)=Wp(Q_succ(i)); 
                end
                Vmax(p)=max(VN);    
             else
                Vmax(p)=2*Hmax;  
            end            
        end
        [v,c]=min(Vmax);
        Wn(q)=min(Vn(q),v); 
        
        if Wn(q)~=Wp(q)
            K(q)=iter;
            C_mF(q,:)=(Vmax==v);
        end
    end             
    fixpoint=isequal(Wp,Wn);
end
disp('Fixed point reached')
Iterations=iter


toc
           
disp('Abstraction and controller synthesis without Nested approach finished')

