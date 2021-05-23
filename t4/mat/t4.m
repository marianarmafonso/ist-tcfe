close all
clear all
pkg load symbolic;

%gain stage

VT=25e-3;
BFN=178.7;
VAFN=69.7;
RE1=100;
RC1=571;
RB1=125050;
RB2=20000;
VBEON=0.7;
VCC=12;
RS=100;
Ci = 4.5e-4;
Cb = 4.465e-3;

%DC analysis
RB=1/(1/RB1+1/RB2);
VEQ=RB2/(RB1+RB2)*VCC;
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1);
IC1=BFN*IB1;
IE1=(1+BFN)*IB1;
VE1=RE1*IE1;
VO1=VCC-RC1*IC1;
VCE=VO1-VE1;



%incremental analysis
gm1=IC1/VT;
rpi1=BFN/gm1;
ro1=VAFN/IC1;
RSB=RB*RS/(RB+RS);

RE1=0;
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2);
AVI_DB = 20*log10(abs(AV1));
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1);
AVIsimple_DB = 20*log10(abs(AV1simple));

ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)));
ZO1 = 1/(1/ro1+1/RC1);

%ouput stage
BFP = 227.3;
VAFP = 37.2;
RE2 = 100;
VEBON = 0.7;
Co = 1.806e-3;

%DC analysis
VI2 = VO1;
IE2 = (VCC-VEBON-VI2)/RE2;
IC2 = BFP/(BFP+1)*IE2;
VO2 = VCC - RE2*IE2;


gm2 = IC2/VT;
go2 = IC2/VAFP;
gpi2 = gm2/BFP;
ge2 = 1/RE2;

AV2 = gm2/(gm2+gpi2+go2+ge2);

ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2);

ZO2 = 1/(gm2+gpi2+go2+ge2);

%total
gB = 1/(1/gpi2+ZO1);
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1;
AV_DB = 20*log10(abs(AV));
ZI=ZI1;
ZO=1/(go2+gm2/gpi2*gB+ge2+gB);

RE1 = 100;
%impedancias para lco (condensadores em curto circuito)
Zx = 1/(1/RE1 + gm1*rpi1/(rpi1+RSB));
Zb = (RE1*Zx)/(RE1+Zx);
lco = (1/(2*pi))*((1/(ZO*Co)) + (1/(ZI*Ci)) + (1/(Zb*Cb)));

n_transistors = 2;
cost1 = (0.1*n_transistors + (Ci + Cb + Co)*10e6 + (RE1 + RE2 + RC1 + RB1 + RB2 + RS)*10e-3);
in = 0;
cost = 68671.21;
in2 = 0;
out = 0;
save("-ascii","../doc/tabelaN.tex", "VEQ", "VO1", "VE1", "VO2", "in", "in2", "out", "VCC");
save("-ascii","../doc/tabelaVal.tex","n_transistors", "Ci", "Cb", "Co", "RE1", "RE2", "RC1", "RB1", "RB2", "RS");
save("-ascii","../doc/tabela2.tex","AV1", "ZI1", "ZO1", "AV2", "ZI2", "ZO2");
save("-ascii","../doc/tabelaRes.tex","cost", "lco", "AV", "ZI", "ZO");
%plot
f = logspace(1, 8, 50);
i = 1;
for (i = 1:50)
 if (lco < f(i))
  gain(i) = AV_DB;
 else
  gain(i) = 20*20*log10(f(i)) + AV_DB - 20*20*log10(lco);
 endif
  i = i + 1;
 endfor
hg = figure();
semilogx (f, gain);
xlabel("f[Hz]");
ylabel("V_o(f)/V_i(f)");
title('Frequency response using logarithmical scale');
print (hg, "../doc/plot1.eps");