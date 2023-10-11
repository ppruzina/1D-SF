function [R0,A0,B0,C0,Q0,e0] = ddcinstabilityjuly22()
syms g d e 
g0 = 1;
guess = 1;
tau = 0.01;
delta = 0.001;
sigma = 10;
epsilon = 1;


R = g/d;
R0 = linspace(1,25,500);
gamma = (R.^2)/7 - R/2 + 5/4;
D = (R*tau-gamma)./(gamma-R);
% D = sqrt(e^2+delta*R^2)/R;

f = D^2*g/(D+1);
c = D^2*d/(D+tau);
% if g0<0
%DC
% p = sigma*(-D^2*g/(D+1) + D^2*d/(D+tau)) - epsilon*e^2/D + mu*(e+delta)/e;

% end
% if g0>0
%SF
p = sigma*(-D^2*g/(D+1) + D^2*d/(D+tau)) - epsilon*e^2/D;
% end
% if g0 == 0
%     disp('choose a non-zero value of $g_0$!')
% end
% 
fg = diff(f,g);
fd = diff(f,d);
fe = diff(f,e);
cg = diff(c,g);
cd = diff(c,d);
ce = diff(c,e);
pg = diff(p,g);
pd = diff(p,d);
pe = diff(p,e);

Fg = (fg*pe-fe*pg)/pe;
Fd = (fd*pe-fe*pd)/pe;
Cg = (cg*pe-ce*pg)/pe;
Cd = (cd*pe-ce*pd)/pe;

A = simplify(Fg*Cd-Fd*Cg);

B = simplify(pe);
C = simplify(Fg+Cd);
Q = simplify(fg*cd-fd*cg);
U = fg*fe*pg + cd*ce*pd + fd*ce*pg + fe*cg*pd;

for i = 1:length(R0)
    p0(i) = subs(p,[g,d],[g0,g0/R0(i)]);
    temp = vpasolve(p0(i),e,guess);
    temp = real(temp(abs(imag(temp))<0.000000001));
    temp = temp(real(temp)>=0);
    if isempty(temp)
        e0(i) = NaN;
    else
        temp = sort(temp);
%         e0(i) = temp(index);  
        e0(i) = max(temp);
        guess = e0(i);
    end
end


l = D/e^(1/2);
for i = 1:length(R0)
    A0(i) = double(subs(A,[g,d,e],[g0,g0/R0(i),e0(i)]));
    B0(i) = double(subs(B,[g,d,e],[g0,g0/R0(i),e0(i)]));
    C0(i) = double(subs(C,[g,d,e],[g0,g0/R0(i),e0(i)]));
    Q0(i) = double(subs(Q,[g,d,e],[g0,g0/R0(i),e0(i)]));
    l0(i) = double(subs(l,[g,d,e],[g0,g0/R0(i),e0(i)]));
    U0(i) = double(subs(U,[g,d,e],[g0,g0/R0(i),e0(i)]));
end
figure()
% a1 = subplot(2,1,1);
plot(R0,A0/abs(median(A0(~isnan(A0)))),'k',R0,-B0/abs(median(B0(~isnan(B0)))),R0,C0/abs(median(C0(~isnan(C0)))),R0,Q0/abs(median(Q0(~isnan(Q0)))),R0,U0/abs(median(U0(~isnan(U0)))),'LineWidth',1)
% plot(R0,A0/abs(median(A0(~isnan(A0)))))
xlabel('$R_0$','interpreter','latex')
hold on
plot([min(R0) max(R0)], [0 0], 'b','LineWidth',1)
hold off
ylim([min(A0/abs(median(A0(~isnan(A0))))),-min(A0/abs(median(A0(~isnan(A0)))))])
% ylim([-4 4])
% legend({'$F_gC_d-F_dC_g$';'$p_e$';'$F_g+C_d$';'$f_gc_d-f_dc_g$'},'interpreter','latex')
% title(['$\epsilon = ',num2str(epsilon),'$, $\delta = ',num2str(delta),'$'],'interpreter','latex')
% a2 = subplot(2,1,2);
% plot(a2,R0,e0,R0,l0)
% xlabel('$R_0$','interpreter','latex')
% legend(a2,{'$e_0$';'$l_0$'},'interpreter','latex')
% title('(b)','interpreter','latex')
% a3 = subplot(2,2,3);
% plot(a3,e0,l0);
% xlabel('$e_0$','interpreter','latex')
% ylabel('$l_0$','interpreter','latex')
% title('(c)','interpreter','latex')
% ylim([-4 4])
set(gcf,'PaperSize',[21 21])
set(gcf,'PaperPosition',[0 0 21 21])
end
