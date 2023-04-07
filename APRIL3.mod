var ct ${C^T}$
    cn ${C^N}$
    n ${N}$
    b ${B}$
    bstar ${B^{*}}$
    bpstar ${B^{p*}}$
    bgstar ${B^{g*}}$
    m ${M}$
    gm ${g_m}$
    gy ${g_y}$
    lambdaa ${\lambda}$
    q ${q}$
    pin ${\pi^N}$
    ptilden ${\tilde{pN}}$
    pii ${\pi}$
    yn ${Y^N}$
    epsiloon ${\varepsilon}$
    ib ${i^b}$
    dmr ${DMR}$
    dmc ${DMC}$
    s ${s}$
    e ${e}$
    psii ${\psi}$
    w ${W}$
;

varexo ua ${u^a}$
;  

parameters ggamma ${\Gamma}$
    deltaa ${\delta}$
    omegaa ${\Omega}$
    psibstar ${\bar{\psi}}$
    thetaa ${\theta}$
    sigmaa ${\sigma}$
    etaa ${\eta}$
    chii ${\chi}$
    betaa ${\beta}$
    alphaa ${\alpha}$
    muu ${\mu}$
    rhog ${\rho_g}$
    ibstar ${i^{b*}}$
    kappaepsilon ${\kappa_{\varepsilon}}$
    kappay ${\kappa_{y}}$
    kappapi ${\kappa_{\pi}}$
    yt ${y^T}$
;

betaa=0.995;
sigmaa=2;
deltaa=0.26;
chii=0.06;
etaa=2;
omegaa=0.6;
psibstar=0.9;
muu=10;
alphaa=0.75;
thetaa=0.7;
ibstar=0.0316;
kappaepsilon = 0.8;
kappay = 0.8;
kappapi = 0.8;
yt=1;
rhog=0.9;
ggamma=deltaa^deltaa*(1-deltaa)^(1-deltaa);
ibstar=0.0316;

model;
ggamma*q^(deltaa-1)*ct+omegaa/(2*epsiloon)*(b(+1)+epsiloon*bpstar(+1))*(psii(+1)-psibstar)^2+bstar(+1)-(1+ibstar)*bstar=ggamma*q^(deltaa-1)*yt;
cn=yn;
thetaa*ct^(deltaa-1-sigmaa*deltaa)*cn^((1-deltaa)*(1-sigmaa))=ggamma*q^(deltaa-1)*lambdaa;
(1-thetaa)*ct^(deltaa-sigmaa*deltaa)*cn^(-sigmaa*(1-deltaa)-thetaa)=ggamma*q^(deltaa)*lambdaa;
n^etaa=lambdaa*w/epsiloon;
chii*(m/epsiloon)^(-1)=lambdaa-betaa*lambdaa(+1)*epsiloon/epsiloon(+1);
betaa*(1+ib(+1))/epsiloon(+1)*lambdaa(+1)=lambdaa*(1/epsiloon*(1+omegaa/2*(psii(+1)-psibstar)^2)+(b(+1)/epsiloon+bpstar(+1))*(omegaa*(psii(+1)-psibstar)*(epsiloon*bpstar)/(b(+1)+epsiloon*bpstar(+1))^2));
betaa*(1+ibstar(+1))*epsiloon(+1)/epsiloon(+1)*lambdaa(+1)=lambdaa*(epsiloon/epsiloon*(1+omegaa/2*(psii(+1)-psibstar)^2)+(b(+1)/epsiloon+bpstar(+1))*(omegaa*(psii(+1)-psibstar)*(-epsiloon*b)/(b(+1)+epsiloon*bpstar(+1))^2));
psii(+1)=b(+1)/(b(+1)+bpstar(+1)*epsiloon);
bstar=bpstar+bgstar;
q=(1-deltaa)/deltaa*ct/cn;
yn=s^(-alphaa)*n^alphaa;
s=thetaa*s(-1)*pin^(muu/alphaa)+(1-thetaa)*ptilden^(-muu/alphaa);
1=thetaa*pin^(muu-1)+(1-thetaa)*ptilden^(1-muu);
dmr=((muu-1)/muu)*q*yn*ptilden^(1-muu)+betaa*thetaa*lambdaa(+1)/lambdaa*(q/q(+1))^(1-deltaa)*(ptilden/ptilden(+1)*1/pin(+1))^(1-muu)*dmr(+1);
dmc=1/alphaa*w*yn^(1/alphaa)*ptilden^(-muu/alphaa)+betaa*thetaa*lambdaa(+1)/lambdaa*(q/q(+1))^(1-deltaa)*(ptilden/ptilden(+1)*1/pin(+1))^(-muu/sigmaa)*dmc(+1);
dmr=dmc;
pii=e^deltaa*pin^(1-deltaa);
epsiloon*(bgstar(+1)-(1+ibstar)*bgstar)-(b(+1)-(1+ib)*b)=m-m(-1);
gm=kappay*gy+kappapi*pii+kappaepsilon*e;
gm=m/m(-1);
gy=yn/yn(-1);
gm=rhog*gm(-1)+ua;
e=epsiloon/epsiloon(-1);
end;

write_latex_static_model;
write_latex_parameter_table;
write_latex_dynamic_model;
collect_latex_files ;

shocks;
var ua;
stderr 0.01;
end;

steady;
%check;
%stoch_simul(order=1,irf=40); 
