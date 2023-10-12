function [smax] = ddcwavenumbergr(R0)
%Produce wavenumber-growth rate plots for the DDC system

syms g d e m s

smax = 0;

tau = 0.0001; % diffusivity ratio
sigma = 10^4; % prandtl number
epsilon = 1; % dissipation parameter
delta = 0.001; % parameter in length scale
W = 0; % optional constant power forcing term
g0 = 1; % background temperature gradient, nondimensionalised to +/- 1
% R0 = 0.8; %density ratio 
Mmax = 0.8; %maximum wavenumber for plot

%set up expressions
R = g/d;
D = sqrt(e^2+delta*R^2)/R; %turbulent diffusivity
f = D^2*g/(D+1); % temperature flux
c = D^2*d/(D+tau); % salinity flux
% generalised source term, with contributions from PE-KE transfer, dissipation, and the constant power source
p = sigma*(-D^2*g/(D+1) + D^2*d/(D+tau)) - epsilon*e^2/D + W; 

%energy diffusivity
k = D^2/(D+sigma) + sigma;

%calculate partial derivatives
fg = diff(f,g);
fd = diff(f,d);
fe = diff(f,e);
cg = diff(c,g);
cd = diff(c,d);
ce = diff(c,e);
pg = diff(p,g);
pd = diff(p,d);
pe = diff(p,e);


%solve steady state energy equation p=0 to find e0(R0)
p0 = subs(p,[g,d],[g0,g0/R0]);
temp = vpasolve(p0,e,0.6574);
if isempty(temp)
    e0 = NaN;
else
    temp = sort(temp);
    temp = temp(abs(imag(temp))<0.000000001);
    temp = temp(temp>=0);
    e0 = max(temp);
end

e0;


%substitute specific values of g0, d0, e0 into partial derivatives
fg0 = double(subs(fg,[g,d,e],[g0,g0/R0,e0]));
fd0 = double(subs(fd,[g,d,e],[g0,g0/R0,e0]));
fe0 = double(subs(fe,[g,d,e],[g0,g0/R0,e0]));
cg0 = double(subs(cg,[g,d,e],[g0,g0/R0,e0]));
cd0 = double(subs(cd,[g,d,e],[g0,g0/R0,e0]));
ce0 = double(subs(ce,[g,d,e],[g0,g0/R0,e0]));
pg0 = double(subs(pg,[g,d,e],[g0,g0/R0,e0]));
pd0 = double(subs(pd,[g,d,e],[g0,g0/R0,e0]));
pe0 = double(subs(pe,[g,d,e],[g0,g0/R0,e0]));
k0 = double(subs(k,[g,d,e],[g0,g0/R0,e0]));


%Linearised system in matrix form Mat = 0
Mat = [s + m^2*fg0, m^2*fd0, m^2*fe0;
       m^2*cg0, s + m^2*cd0, m^2*ce0;
       -pg0, -pd0, s + m^2*k0 - pe0];
%find characteristic polynomial as determinant of Mat
spoly = det(Mat);

%set up array of wavenumbers and blank array for solutions
M = linspace(0.0,Mmax,500);
ssol = zeros(3,length(M));

%solve characteristic polynomial for growth rate s at each value of M
for i = 1:length(M)
    spolym = subs(spoly,m,M(i));
    temp = solve(spolym,s);
    if ~isempty(temp)
        ssol(:,i) = temp;
    else
        ssol(:,i) = [NaN;NaN;NaN];
    end
end
%filter for only positive, real values of s
ssol = sort(ssol,1,'ComparisonMethod','real');
spos = ssol(ssol>=0);

if ~isempty(spos)
    sposi = ssol;
    sposi(sposi<0) = NaN;
    smaxes = max(sposi,[],1,'omitnan');
else
smaxes = max(ssol,[],1);
end

%find maximum growth rate smax, and corresponding wavenumber mmax
[smax,i] = max(smaxes,[],'All','ComparisonMethod','real');
mmax = M(i)





% Create figure showing the growth rates as functions of wavenumber.

% F = figure();
% subplot(1,2,1);
plot(M,ssol,'k')
xlabel('wavenumber $m$','interpreter','latex')
ylabel('growth rate $s$','interpreter','latex')
if smax>0
ylim([-2*real(smax) 1.2*real(smax)])
end
yline(0);
% title('(d)','interpreter','latex')
% subplot(1,2,2)
% plot(M,imag(ssol),'k')
end