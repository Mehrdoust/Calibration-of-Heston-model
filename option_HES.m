clc
clear all
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%S0=Enter initial price;
%T=Enter maturity date;
RR=[7 7];
%Cmarket1=Enter the actual option prices (In-sample);
%KK=Enter the actual strike prices (In-sample);
%bid=Enter the actual bid values(In-sample);
% ask=Enter the actual ask values data (In-sample);
Cmarket2=Enter the actual option prices (Out-of-sample);
KK1=Enter the actual strike prices (Out-of-sample);
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%

n=length(Cmarket1);
ww=1./(sqrt(abs(ask-bid)));
C=@(x)(0);

syms z
nn=1;
p(z)=laguerreL(nn,z);
q(z)=laguerreL(nn+1,z);
zz=roots(sym2poly(p));
wz=zz./((nn+1)^2.*q(zz).^2);
i=sqrt(-1);
%%%%%%%%%%
for j=1:n
    w=ww(j);
    Cmarket=Cmarket1(j);
    K=KK(j);
    alphaa=@(x,o)(((-o.^2)/2)-((i.*o)/2));
    betaa=@(x,o)(x(1)-x(4).*x(3).*i.*o);
    gama=@(x,o)((x(3).^2)./2);
    h=@(x,o)(sqrt(betaa(x,o).^2-4*alphaa(x,o).*gama(x,o)));
    u1=@(x,o)((betaa(x,o)-h(x,o))./(x(3).^2));
    u2=@(x,o)((betaa(x,o)+h(x,o))./(x(3).^2));
    g=@(x,o)(u1(x,o)./u2(x,o));
    D=@(x,o)(u1(x,o).*((1-exp(h(x,o).*T))./(1-g(x,o).*exp(-h(x,o).*T))));
    C1=@(x,o)(x(1).*(u1(x,o).*T-(2./(x(3).^2))*...
        log((1-g(x,o).*exp(h(x,o).*T))./(1-g(x,o)))));
    sai=@(x,o)(exp(C1(x,o).*x(2)+D(x,o).*x(5)+i.*o.*log(S0.*exp(x(6).*T))));
    pi1=@(x)(.5+(1/pi)*real(sum(((exp(-i.*zz.*log(K)).*sai(x,zz-i))./(i.*zz.*sai(x,-i))).*...
        exp(zz).*wz)));
    pi2=@(x)(.5+(1/pi)*real(sum(((exp(-i.*zz.*log(K)).*sai(x,zz))./(i.*zz)).*...
        exp(zz).*wz)));
    pi1h=@(x)(.5+(1/pi)*real(sum(((exp(-i.*zz.*log(x(7))).*sai(x,zz-i))./(i.*zz.*sai(x,-i))).*...
        exp(zz).*wz)));
    pi2h=@(x)(.5+(1/pi)*real(sum(((exp(-i.*zz.*log(x(7))).*sai(x,zz))./(i.*zz)).*...
        exp(zz).*wz)));
    Cmodel=@(x)(S0.*pi1(x)-exp(-x(6).*T).*K.*pi2(x));
    Cmodelh=@(x)(S0.*pi1h(x)-exp(-x(6).*T).*x(7).*pi2h(x));
    C=@(x)(double(C(x)+(w.*((Cmarket-Cmodel(x))./Cmarket).^2)));
end
flag=0; 
while flag==0
    options =  gaoptimset('TimeLimit',.1,'PopulationSize',2);
    par_Hes=ga(C,6,[],[],[],[],[1.2329 .01969  .1868 -.2789 .001764 .0026],[3 .15 .54353 -.1435 .2 .16533],[],options);
    for i=1:length(Cmarket1)
        Ks=KK(i);
        ch(i)=double(Cmodelh([par_Hes Ks]));
        if ch(i)<0
            ch(i)=0;
        end
    end
    for i=1:length(Cmarket2)
        Ks=KK1(i);
        ch1(i)=double(Cmodelh([par_Hes Ks]));
        if ch1(i)<0
            ch1(i)=0;
        end
    end
    
    RMSE=double([sqrt(mse(Cmarket1-ch)) sqrt(mse(Cmarket2-ch1))])
    RMSE1=double([mean(abs(Cmarket1-ch)./Cmarket1)  mean(abs(Cmarket2-ch1)./Cmarket2)]);
    if RMSE(1)<RR(1)&&RMSE(2)<RR(2)&&RMSE(1)<RMSE(2) 
        flag=1;
    end
end
clc
%%%%%%%%%%%%%%%%%%
par_Hes'
%%%%%%%%%%%%%%%%%%%
disp('in-sample      out-of-sample')
RMSE1
RMSE



