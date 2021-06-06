close all
clear all
pkg load symbolic;

C1 = 220e-9;
C2 = 220e-9;
R1 = 1e3;
R2 = 1e3;
R3 = 100e3;
R4 = 1e3;

f = logspace(1,8,70);
w = 2*pi*f;
Zc1 = 1./(j*w*C1);
Zc2 = 1./(j*w*C2);
%transfer functions
Tf = (R1*C1*w*i)./(1+R1*C1*w*i)*(1+R3/R4).*(1./(1+R2*C2*w*i));

%frequency response
hf = figure();
semilogx(f,20*log10(abs(Tf)));
xlabel("Frequency [Hz]");
ylabel("Gain [dB]");
title("Gain");
print(hf, "../doc/gain.eps", "-depsc");

hg = figure();
semilogx(f,180*arg(Tf)/pi);
xlabel("Frequency [Hz]");
ylabel("Phase [Degrees]");
title("Phase");
print(hg, "../doc/phase.eps", "-depsc");

%central frequency
wL = 1/(R1*C1);
lf = wL/(2*pi);
wH = 1/(R2*C2);
fh = wH/(2*pi);
wO = sqrt(wL*wH);
cf = wO/(2*pi);
cf2 = 1000;
w2 = 2*pi*cf2;

Gain = abs((R1*C1*w2*i)/(1+R1*C1*w2*i)*(1+R3/R4)*(1/(1+R2*C2*w2*i)));
Gain_db = 20*log10(abs(Gain));

%output e input impedances na central frequency
ZC1 = 1/(i*w2*C1);
ZC2 = 1/(i*w2*C2);
ZI = abs(ZC1 + R1); %C1 e R1 em serie
ZO = abs((R2*ZC2)/(R2+ZC2)); %C2 e R2 em paralelo

save("-ascii","../doc/tabelaVal.tex","R1", "R2", "R3", "R4", "C1", "C2");
save("-ascii","../doc/tabelaRes.tex","lf", "fh", "cf", "ZI", "ZO", "Gain", "Gain_db");
