clear;
clc;
rng(0)
%{
re表示MC仿真重复次数
T代表模拟样本的长度
dt代表去趋势方式：'1'代表一步去趋势；"2"代表两步去趋势，"3"代表GLS去趋势
q表示去趋势阶数：q=0表示去常数；q=1表示去线性趋势
c表示偏离单位根程度：a=1-c/T；
var表示构造MC样本时扰动项方差；
B表示使用bootstrap方法时构造的bootstrap样本数量；
delta表示构造MC样本时趋势项的系数向量；
pp表示Z检验的方法：“1”表示传统的PP检验，即使用Newey-West进行长期方差估计；“2”表示经过修正的PP单位根检验，即使用固定平滑参数的HAR长期方差估计代替Newey-West长期方差估计，同时对统计量进行一定的修正，使得统计量极限分布与长期方差无关；“3”表示使用Wild Bootstrap方法获得修正PP检验的临界值；“4”表示使用Sieve Bootstrap方法获得修正PP检验的临界值；
b 表示固定平滑参数的HAR长期方差估计的截断参数：在“type”为“bt”、“pr”、“qs” 时，b表示截断参数；在 “type”为“os”时表示使用的正交序列的长度；
type 表示估计长期方差协方差矩阵时使用的核函数种类或非参数方法：“bt”表示“Bartlett kernel”；“pr”表示“Parzen kernel”；“qs”表示“QuadraticSpectral kernel”；“os”表示“Orthonormal Series”方法；
%}

re=2000;

c=-5;
delta=[0;0];
var=1;
B=499;
ic='aic';

b=0.02;
%num_bs=length(bs);

qs=[-1 0 1];

T=100;

Bs=[99 199 299 399 499];
rhos=[0, 0.5, -0.3, 0, 0;
    0,  0,  -0.2, 0, 0];
phis=[0, 0, 0, -0.5, 0.3;
    0, 0, 0,   0,  0.2];
pp=4;
dt=3;
type='bt';
size_zas=zeros(5,12);
size_zts=zeros(5,12);

    for e=1:3 % Nodetrended or Demeaned or Detrend
        q=qs(e);
        for f=1:5  % Constant or AR(1) or MA(1)
            rho=rhos(:,f)';
            phi=phis(:,f)';
            for g=1:5 % PP with OLS_detrend or PP(fb) with OLS detrend or PP^b(fb) with GLS detrend
                B=Bs(g);
                
                zalpha=zeros(re,1);
                p_zalpha=zeros(re,1);
                zt=zeros(re,1);
                p_zt=zeros(re,1);
                
                
                y=DGP(T,re,c,delta,rho,phi,var);
                parfor j=1:re
                    s = RandStream('mt19937ar','Seed',j);
                    RandStream.setGlobalStream(s);
                    [zalpha(j,1),p_zalpha(j,1),zt(j,1),p_zt(j,1)]=ztest(y(:,j),dt,q,pp,ic,B,b,type);
                end
                
                size_za=sum(p_zalpha(:,1)<=0.05)/re;
                size_zt=sum(p_zt(:,1)<=0.05)/re;
                size_zas(f,(e-1)*5+g)=size_za;
                size_zts(f,(e-1)*5+g)=size_zt;
            end
        end
    end
save('PowerTableforT_100withDifferentBootsam','size_zas','size_zts');
