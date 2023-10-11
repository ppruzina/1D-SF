function [smax] = ddcwavenumbergr(R0)
syms g d e m s

smax = 0;

tau = 0.0001;
sigma = 10^4;
epsilon = 1;
delta = 0.001;
W = 0;
g0 = 1;
% R0 = 0.8;
Mmax = 0.8;
R = g/d;
D = sqrt(e^2+delta*R^2)/R;
f = D^2*g/(D+1);
c = D^2*d/(D+tau);

p = sigma*(-D^2*g/(D+1) + D^2*d/(D+tau)) - epsilon*e^2/D + W;


k = D^2/(D+sigma) + sigma;


fg = diff(f,g);
fd = diff(f,d);
fe = diff(f,e);
cg = diff(c,g);
cd = diff(c,d);
ce = diff(c,e);
pg = diff(p,g);
pd = diff(p,d);
pe = diff(p,e);



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

%matrix to find growth rate equation
Mat = [s + m^2*fg0, m^2*fd0, m^2*fe0;
       m^2*cg0, s + m^2*cd0, m^2*ce0;
       -pg0, -pd0, s + m^2*k0 - pe0];
spoly = det(Mat);

M = linspace(0.0,Mmax,500);

ssol = zeros(3,length(M));
for i = 1:length(M)
    spolym = subs(spoly,m,M(i));
    temp = solve(spolym,s);
    if ~isempty(temp)
        ssol(:,i) = temp;
    else
        ssol(:,i) = [NaN;NaN;NaN];
    end
end
ssol = sort(ssol,1,'ComparisonMethod','real');
spos = ssol(ssol>=0);

if ~isempty(spos)
    sposi = ssol;
    sposi(sposi<0) = NaN;
    smaxes = max(sposi,[],1,'omitnan');
else
smaxes = max(ssol,[],1);
end
[smax,i] = max(smaxes,[],'All','ComparisonMethod','real');
mmax = M(i)






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