function [sol,xmesh,tmesh] = ddcsolverjuly22(tend)
%Function solves the Salt fingering model
tic
%Inputs: 
% n: initial number of sinusoidal perturbation
% R0: initial density ratio
%s: = +/- 1, sets Tz in SF or DC
% solold: previous solution to be used as i.c. Set as 1 to use uniform gradient plus sine perturbation
% tstart: starting time. If using solold, take previous tend value
% tend: end time of integration
% Peinv: inverse Peclet number, controls molecular diffusion
% Reinv: inverse Reynolds number, controls viscous diffusion
% H: domain depth
% epsilon: dissipation parameter

timeseries = 0;%Set to 1 to plot a timeseries of bz
endplot = 1; %Set to 1 to plot the final state
growthcheck = 1;
solold = 1;

%Parameters can be set here or as inputs to the function
%Physical parameters

tau = 0.0001;
sigma = 10^4;
epsilon = 1;
delta = 0.001;
W = 0;

R0 = 3.64;
g0 = 1;
d0 = g0/R0;
%Domain size and mesh parameters
H = 2000;
n = 6;
% tend = 1000000;
Nx = 4*H;
% Nt =10000;


evar = sym('evar');
D0 =  sqrt(evar^2+delta*R0^2)/R0;


p = -sigma*(D0^2*g0/(D0+1) - D0^2*d0/(D0+tau)) - epsilon*evar^2/D0 + W; 


e0 = max(vpasolve(p,evar,0.2))
Hyperdiffusion = 0;%Set as 1 to add hyperdifussion to the system




gpert = g0/1000; %perturbation gradient
m = 2*pi*n/H

%Find coefficients for the perturbations.
%Mat is the linear matrix according to the RHS of the equations. find
%eigenvalues D, corresponding to growth rates, then we will initialise with
%perturbations in T,S and e with coefficients according to eigenvector of
%maximum growth


Mat = [-(m^2*(delta*g0^2 - d0^2*e0^2 + d0*delta*g0*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2)))/(delta*g0^2 + g0^2 + d0^2*e0^2 + 2*d0*g0*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2)), -(m^2*(d0^2*((d0^2*e0^2 + delta*g0^2)/d0^2)^(3/2) - delta*g0^2*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2) + 2*d0*e0^2*g0))/(delta*g0^2 + g0^2 + d0^2*e0^2 + 2*d0*g0*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2)), -(d0^2*e0*m^2*(2*g0 + d0*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2)))/(delta*g0^2 + g0^2 + d0^2*e0^2 + 2*d0*g0*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2))
(d0^2*m^2*(d0^2*((d0^2*e0^2 + delta*g0^2)/d0^2)^(3/2) - delta*g0^2*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2) + 2*d0*e0^2*g0*tau))/(g0^2*(delta*g0^2 + d0^2*e0^2 + g0^2*tau^2 + 2*d0*g0*tau*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2))), -(m^2*(2*d0^3*((d0^2*e0^2 + delta*g0^2)/d0^2)^(3/2) + delta*g0^3*tau - d0*delta*g0^2*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2) + 3*d0^2*e0^2*g0*tau))/(delta*g0^3 + g0^3*tau^2 + d0^2*e0^2*g0 + 2*d0*g0^2*tau*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2)), -(d0^3*e0*m^2*(d0*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2) + 2*g0*tau))/(delta*g0^3 + g0^3*tau^2 + d0^2*e0^2*g0 + 2*d0*g0^2*tau*((d0^2*e0^2 + delta*g0^2)/d0^2)^(1/2))
(delta*e0^2*epsilon*g0^2)/(d0^3*(e0^2 + (delta*g0^2)/d0^2)^(3/2)) - (e0^2*epsilon)/(d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2)) - sigma*((2*delta)/((d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0 + 1) - (d0^2*(e0^2 + (delta*g0^2)/d0^2))/(g0^2*((d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0 + 1)) - (2*d0*delta)/(g0*(tau + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0)) + (2*d0^3*(e0^2 + (delta*g0^2)/d0^2))/(g0^3*(tau + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0)) - (d0^2*(delta/(d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2)) - (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0^2)*(e0^2 + (delta*g0^2)/d0^2))/(g0*((d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0 + 1)^2) + (d0^3*(delta/(d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2)) - (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0^2)*(e0^2 + (delta*g0^2)/d0^2))/(g0^2*(tau + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0)^2)), (e0^2*epsilon*g0)/(d0^2*(e0^2 + (delta*g0^2)/d0^2)^(1/2)) - sigma*((2*delta)/(tau + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0) - (2*delta*g0)/(d0*((d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0 + 1)) + (2*d0*(e0^2 + (delta*g0^2)/d0^2))/(g0*((d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0 + 1)) - (3*d0^2*(e0^2 + (delta*g0^2)/d0^2))/(g0^2*(tau + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0)) - (d0^2*((e0^2 + (delta*g0^2)/d0^2)^(1/2)/g0 - (delta*g0)/(d0^2*(e0^2 + (delta*g0^2)/d0^2)^(1/2)))*(e0^2 + (delta*g0^2)/d0^2))/(g0*((d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0 + 1)^2) + (d0^3*((e0^2 + (delta*g0^2)/d0^2)^(1/2)/g0 - (delta*g0)/(d0^2*(e0^2 + (delta*g0^2)/d0^2)^(1/2)))*(e0^2 + (delta*g0^2)/d0^2))/(g0^2*(tau + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0)^2)) - (delta*e0^2*epsilon*g0^3)/(d0^4*(e0^2 + (delta*g0^2)/d0^2)^(3/2)), (e0^3*epsilon*g0)/(d0*(e0^2 + (delta*g0^2)/d0^2)^(3/2)) - m^2*(sigma + (d0^2*(e0^2 + (delta*g0^2)/d0^2))/(g0^2*(sigma + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0))) - (2*e0*epsilon*g0)/(d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2)) - sigma*((2*d0^2*e0)/(g0*((d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0 + 1)) - (2*d0^3*e0)/(g0^2*(tau + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0)) - (d0^3*e0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/(g0^2*((d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0 + 1)^2) + (d0^4*e0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/(g0^3*(tau + (d0*(e0^2 + (delta*g0^2)/d0^2)^(1/2))/g0)^2))];

 

[V,D] = eig(Mat);
D = diag(D);
D(abs(imag(D))>0.0000001) = NaN;
[spredict,Igr] = max(D);
spredict = double(spredict)
dpert = real(gpert*V(2,Igr)/V(1,Igr));
epert = real(gpert*V(3,Igr)/V(1,Igr));




global counter
counter = 0; 

tstart = 0;
xmesh = linspace(0,H,Nx);
if solold == 1
    tmesh = logspace(0,log10(tend),log10(tend)*100+1);
    tmesh = [0 tmesh(3:end)];
else
tmesh = logspace(log10(tstart),log10(tend),(log10(tend)-log10(tstart))*500);
end
%Set tolerance
tol = 1e-5  ;
options = odeset('RelTol',tol,'AbsTol',1e-3*tol);
sol = pdepe(0,@pde,@(x)ics(x,solold),@bcs,xmesh,tmesh,options);
% end

Esol = sol(:,:,3);
Tsol = sol(:,:,1);
Ssol = sol(:,:,2);
Ntn = size(sol,1);
%calculate gradient field
Tzsol = diff(Tsol,1,2)./diff(xmesh);
Tzsol = [Tzsol(:,1) Tzsol];%repeat first entry to keep length of Nx
Szsol = diff(Ssol,1,2)./diff(xmesh);
Szsol = [Szsol(:,1) Szsol];




%Plot just the buoyancy gradient (and/or energy) at the start and end 
if endplot == 1
figure(1)
plot(xmesh,Esol(1,:),'--');
% plot(xmesh,10*(Tzsol(1,:)-Szsol(1,:)),'--')
hold on
% set(gca,'ColorOrderIndex',1)
plot(xmesh,Esol(end,:));
% plot(xmesh,10*(Tzsol(end,:)-Szsol(end,:)));
xlabel('$z$','Interpreter','Latex')
title(['$t = $',num2str(tend),', $g_i = $',num2str(g0),', $d_i = $',num2str(d0)],'Interpreter','Latex');
set(gcf,'PaperSize',[15 15])
set(gcf,'PaperPosition',[0 0 15 15]);

end

%Plot a timeseries of the gradient at ten times
if timeseries == 1
B(1,:) = sol(1,:,1) - sol(1,:,2);
for i = 1:10
B(i+1,:) = sol(Nt*i/10,:,1) - sol(Nt*i/10,:,2);
end
for i = 1:11
Bz(i,:) = diff(B(i,:))./diff(xmesh);
end
for i = 1:11
maxes(i) = max(Bz(i,:));
mins(i) = min(Bz(i,:));
Bz(i,:) = (Bz(i,:)-mins(i))/(maxes(i)-mins(i))-1/2;
% end
% cummax = cumsum(maxes);
Bzplus(i,:) = Bz(i,:) + i;
end
xminus = 0.5*(xmesh(1:end-1)+xmesh(2:end));
figure()
plot(xminus,Bzplus,'b')
xlabel('$z$','Interpreter','Latex')
ylabel('$b_z$','Interpreter','Latex')
set(gcf,'PaperSize',[15 15])
set(gcf,'PaperPosition',[0 0 15 15]);
end

if growthcheck ==1
    Eend = Esol(end,:);
    Estart = Esol(1,:);
    Eprev = Esol(end-1,:);
    
    endamp = max(Eend)-min(Eend);
    startamp = max(Estart)-min(Estart);
    prevamp = max(Eprev)-min(Eprev);
    growthrate = log(endamp/startamp)/(tmesh(end))
    localgrowthrate = log(endamp/prevamp)/(tmesh(end)-tmesh(end-1))
    predictgrowthrate = spredict
end

%code pde function
function [c, f, s] = pde(x,t,u,dudx)
    T = u(1);
    S = u(2);
    e = u(3);
    Tz = dudx(1);
    Sz = dudx(2);
    ez = dudx(3);
    R = Tz/Sz;




    %lengthscale

    D = sqrt(e^2 + delta*R^2)/R;
    % %with diffusion
        c = [1;1;1];
        f = [D^2*Tz/(D+1); D^2*Sz/(D+tau); D^2*ez/(D+sigma) + sigma*ez];
        s = [0;0;sigma*(-D^2*Tz/(D+1) + D^2*Sz/(D+tau)) - epsilon*e^2/D + W]; 

    end
 
%coe i.c.s
function u0 = ics(x,solold)
    if solold == 1%if solold = 1, start from initial condition
        %BLY
%     u0 = [(x-H*sinh(20*(x/H - 1/2))/(20*cosh(10)));(1/R0)*(x-H*sinh(20*(x/H - 1/2))/(20*cosh(10)));e0];

%i.c. with m sine modes perturbing uniform gradient

    finit = sin(n*2*pi*x/H);
    finitdiv = 2*pi*n*cos(n*2*pi*x/H)/H;%to initialise with energy
    u0 = [g0*x + gpert*finit;g0*x/R0 + dpert*finit;e0 + epert*finitdiv];



    else%if solold is a previous solution of size(5,Nx), use this as initial condition 
        counter = counter + 1;
        u01 = solold(counter,1);
        u02 = solold(counter,2);
        u03 = solold(counter,3);

       u0 = [u01;u02;u03];
    end



end
%code b.c.s
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)
    %BLY
%         pl = [0;0;0];
%         ql = [1;1;1];
%         pr = [0;0;0];
%         qr = [1;1;1];

%     pl = [ul(1);ul(2);0];
%     ql = [0;0;1];
%     pr = [ur(1) - H;ur(2) - H/R0;0];
%     qr = [0;0;1];

    pl = [ul(1);ul(2);ul(3)-e0];
    ql = [0;0;0];
    pr = [ur(1) - H*g0;ur(2) - H*g0/R0;ur(3)-e0];
    qr = [0;0;0];

    
end




toc
end

