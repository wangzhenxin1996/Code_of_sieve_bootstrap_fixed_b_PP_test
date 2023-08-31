clear
clc

rng(0)

data=readtable('USINFQP.xls');
inflation=data.CPALTT01USQ657N;

dt=3;pp=4;
q=-1;
ic='aic';
B=499;
b=0.02;
type='bt';
[zalpha,p_zalpha,zt,p_zt]=ztest(inflation,dt,q,pp,ic,B,b,type);
dtinf=detrend(inflation,0,3,0);
r=dtinf;
model=arima(1,1,1);
model.Constant=0;
c=estimate(model,r);