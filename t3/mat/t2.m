close all
clear all
pkg load symbolic
pkg load control

values = dlmread('data.txt');
R1 = values(2,4) * (10^3);
R2 = values(3,3) * (10^3);
R3 = values(4,3) * (10^3);
R4 = values(5,3) * (10^3);
R5 = values(6,3) * (10^3);
R6 = values(7,3) * (10^3);
R7 = values(8,3) * (10^3);
Vs = values(9,3);
C = values(10,3) * (10^-6);
Kb = values(11,3) * (10^-3);
Kd = values(12,3) * (10^3);
save("-ascii","../doc/tabelaVal.tex","R1", "R2", "R3", "R4", "R5", "R6", "R7", "Vs", "C", "Kb", "Kd");
R= [R1, R2, R3, R4, R5, R6, R7];
for k=1: 7
    G(k)=1/R(k);
    k=k+1;
end
Id = 0.00102475824097;
Va = 5.13988034104;
Kc = 8161.13797582;
Z=0;
O=1;

% Método dos nós (ponto1)
A= [-G(1), G(1) + G(2) + G(3), -G(2), -G(3), Z, Z, Z;
    Z, -G(2) - Kb, G(2), Kb, Z, Z, Z;
    Z, Kb, Z, -Kb - G(5), G(5), Z, Z;
    Z, Z, Z, Z, Z, G(6)+G(7), -G(7);
    O, Z, Z, Z, Z, Z, Z;
    Z, Z, Z, O, Z, Kd * G(6), -O;
    G(1), -G(1), Z, -G(4), Z, -G(6), Z];
B= [Z; Z; Z; Z; Vs; Z; Z];
V=A\B;
V1 = V(1); 
V2 = V(2); 
V3 = V(3); 
V5 = V(4); 
V6 = V(5); 
V7 = V(6); 
V8 = V(7);
Ib = Kb * (V2 - V5);
IR1 = (V1 - V2) * G(1);
IR2 = (V2 - V3) * G(2);
IR3 = (V2 - V5) * G(3);
IR4 = (-V5) * G(4);
IR5 = (V5 - V6) * G(5);
IR6 = (-V7) * G(6);
IR7 = (V7 - V8) * G(7);
IVs = -IR1;
IVd = IR3 - IR5 + IR4;

save("-ascii","../doc/tabela1.tex","V1", "V2", "V3", "V5", "V6", "V7", "V8", "Ib", "IR1", "IR2", "IR3", "IR4", "IR5", "IR6", "IR7", "IVs", "IVd");

% Ponto2
Vx =  V6 - V8;
E = [G(1) + G(2) + G(3), -G(2), -G(3), Z, Z, Z;
    -G(2) - Kb, G(2), Kb, Z, Z, Z;
    Z, Z, Z, Z, G(6)+G(7), -G(7);
    Z, Z, O, Z, Kd * G(6), -O;
    -G(1), Z, -G(4), Z, -G(6), Z;
    Z, Z, Z, O, Z, -O];
F = [Z; Z; Z; Z; Z; Vx];
Vt0 = E\F;
V2t0 = -Vt0(1);
V3t0 = -Vt0(2);
V5t0 = Vt0(3);
V6t0 = Vt0(4);
V7t0 = Vt0(5);
V8t0 = -Vt0(6);
Ibt0 = Kb * (V2t0 - V5t0);
IR1t0 = ( V2t0) * G(1);
IR2t0 = (V2t0 - V3t0) * G(2);
IR3t0 = (V2t0 - V5t0) * G(3);
IR4t0 = (V5t0) * G(4);
IR5t0 = (V5t0 - V6t0) * G(5);
IR6t0 = (V7t0) * G(6);
IR7t0 = (V7t0 - V8t0) * G(7);
IVst0 = IR1t0;
IVdt0 = IR3t0 - IR5t0 + IR4t0;

Ix = -(V6t0 - V5t0)*G(5) - Kb*(V2t0 - V5t0);
Req = abs(Vx/Ix);
tau = Req * C;
save("-ascii","../doc/tabela2.tex","V2t0","V3t0","V5t0","V6t0","V7t0","V8t0","Vx", "Ix", "Req", "tau", "Ibt0",
"IR1t0","IR2t0", "IR3t0", "IR4t0","IR5t0", "IR6t0", "IR7t0", "IVst0", "IVdt0");

% ponto 3
syms t;
syms v6_n(t);
t=0:1e-6:20e-3;
v6_n = (V8t0 + Vx)*(exp((-t)/tau));
%time axis: 0 to 20ms with 1us steps
hf = figure ();
plot (t*1000, v6_n, "b");
xlabel ("t[ms]");
ylabel ("v_{6n}(t) [V]");
print (hf, "../doc/v6n.eps");
%print v6n.eps

% ponto 4
f= 1000; %Hz
w= 2*(pi)*f;
Zc= 1/(i*w*C);
Yc= 1/Zc;

u = 1;
X= [-G(1), G(1) + G(2) + G(3), -G(2), -G(3), Z, Z, Z;
    Z, -G(2) - Kb, G(2), Kb, Z, Z, Z;
    Z, Kb, Z, -Kb - G(5), G(5) + Yc , Z, -Yc;
    Z, Z, Z, Z, Z, G(6)+G(7), -G(7);
    O, Z, Z, Z, Z, Z, Z;
    Z, Z, Z, O, Z, Kd * G(6), -O;
    G(1), -G(1), Z, -G(4), Z, -G(6), Z];
Y = [Z; Z; Z; Z; u; Z; Z];
V_amplitude = X\Y; % solução com phasor

for n=1:7
      V_magnitude(n)= abs(V_amplitude(n));
      variavel= imag(V_amplitude(n))/real(V_amplitude(n));
      V_phase(n)= atan(variavel);
      n=n+1;
end
V1_magnitude = V_magnitude(1);
V2_magnitude = V_magnitude(2);
V3_magnitude = V_magnitude(3);
V5_magnitude = V_magnitude(4);
V6_magnitude = V_magnitude(5);
V7_magnitude = V_magnitude(6);
V8_magnitude = V_magnitude(7);
V1_phase = V_phase(1);
V2_phase = V_phase(2);
V3_phase = V_phase(3);
V5_phase = V_phase(4);
V6_phase = V_phase(5);
V7_phase = V_phase(6);
V8_phase = V_phase(7);
save("-ascii","../doc/tabela6.tex","V1_magnitude", "V2_magnitude", "V3_magnitude", 
"V5_magnitude", "V6_magnitude", "V7_magnitude", "V8_magnitude", "V1_phase", "V2_phase", "V3_phase", 
"V5_phase", "V6_phase", "V7_phase", "V8_phase");

syms vs(t);
syms Vf6(t);

t=0:1e-6:20e-3;
%phase= atan(w*Req*C) - (pi)/2 ;
Vf6 = V_amplitude(5) * sin(w*t); %solução forçada
t=0:1e-6:20e-3;
vs = u * sin(w*t);

% ponto 5
for(m=1:25001)
      if m < 5001
        v6(m) = V6;
        v_s(m) = Vs;
       elseif m == 5001
        v6(m) = V6t0;
        v_s(m) = Vs;
       elseif m > 5001
        v6(m) = v6_n(m-5001) + (Vf6(m-5001));
        v_s(m) = vs(m-5001);
       end
end
t=-5e-3:1e-6:20e-3;
hg =  figure();
plot (t*1000, v6, "b");
hold on;
plot (t*1000, v_s, "g");
legend("v6","v_s");
xlabel('t[ms]');
ylabel('V');
title('Voltages as functions of time: v_s and v_6 (Volts)');
print (hg, "../doc/v6.eps");

% ponto 6
phasor_vs = 1*power(e,-i*pi/2);
f6 = -1:1e-03:6; %Hz
W = 2*pi*power(10,f6);
Y_c = i .* W .* C;
Z_c = 1. ./Y_c;

X6 = [-Kb - G(2), G(2), Kb, Z;         
     G(3)-Kb,  Z, Kb-G(3)-G(4), -G(6);
     Kb-G(1)-G(3), Z, G(3)-Kb, Z;
     Z, Z, O, Kd * G(6)- R7 * G(6)-O]; 
y6 = [Z; Z; -phasor_vs * G(1); Z];

V6_maluco=X6\y6; % v2, v3, v5, v7
 
v8 = R7*(G(1)+G(6))*V6_maluco(4) + Z*Z_c;
v66 = ((G(5)+Kb)*V6_maluco(3)-Kb*V6_maluco(1)+ (v8 ./ Z_c)) ./ (G(5) + 1. ./ Z_c);
vc6 = v66 - v8;
vs6 = power(e,i*pi/2) + Z*W;

phase_v66 = 180/pi*(angle(v66)); %in degrees
for  var=1:length(phase_v66)
	if(phase_v66(var)<= -90) 
		phase_v66(var) = phase_v66(var) + 180;
	elseif (phase_v66(var)>= 90) 
		phase_v66(var) = phase_v66(var) - 180;
endif
endfor


hg = figure ();
plot (f6, phase_v66, "b");
hold on;
plot (f6, 180/pi*(angle(vc6) + pi), "g");
hold on;
plot (f6, (180*angle(vs6))/pi, "m");

legend("v_6","v_c","v_s");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Phase v_c(f), v_6(f), v_s(f) [degrees]");
print (hg, "../doc/Phase(degrees).eps");

hj = figure ();
plot (f6, 20*log10(abs(v66)), "b");
hold on;
plot (f6, 20*log10(abs(vs6)), "g");
hold on;
plot (f6, 20*log10(abs(vc6)), "m");

legend("v_6","v_s","v_c");
xlabel ("log_{10}(f) [Hz]");
ylabel ("Magnitude v_c, v_6, v_s [dB]");
print (hj, "../doc/MagnitudedB.eps");

%ngspice

filename = 'ngspice_t21.txt'
file = fopen(filename, 'w');
fprintf (file, "Vs 1 0 DC %.11e\nR1 1 2 %.11e\nR2 2 3 %.11e\nR3 2 5 %.11e\nR4 0 5 %.11e\nR5 5 6 %.11e\nR6 9 7 %.11e\nR7 7 8 %.11e\nVaux 0 9 0V\nHd 5 8 Vaux %.11e\nGb 6 3 (2,5) %.11e\nC1 6 8 %.11e", Vs, R1, R2, R3, R4, R5, R6, R7, Kd, Kb, C);
fflush(filename);
fclose(filename);

filename = 'ngspice_t22.txt'
file = fopen(filename, 'w');
fprintf (file, "Vs 1 0 DC 0\nR1 1 2 %.11e\nR2 2 3 %.11e\nR3 2 5 %.11e\nR4 0 5 %.11e\nR5 5 6 %.11e\nR6 9 7 %.11e\nR7 7 8 %.11e\nVaux 0 9 0V\nHd 5 8 Vaux %.11e\nGb 6 3 (2,5) %.11e\nVx 6 8 DC %.11e", R1, R2, R3, R4, R5, R6, R7, Kd, Kb, Vx);
fflush(filename);
fclose(filename);

filename = 'ngspice_t23.txt'
file = fopen(filename, 'w');
fprintf (file, "Vs 1 0 DC 0\nR1 1 2 %.11e\nR2 2 3 %.11e\nR3 2 5 %.11e\nR4 0 5 %.11e\nR5 5 6 %.11e\nR6 9 7 %.11e\nR7 7 8 %.11e\nVaux 0 9 0V\nHd 5 8 Vaux %.11e\nGb 6 3 (2,5) %.11e\nC1 6 8 %.11e\n\n.op\n\n.ic v(6) = %.11e v(8) = 0\n\n.end", R1, R2, R3, R4, R5, R6, R7, Kd, Kb, C, Vx);
fflush(filename);
fclose(filename);

filename = 'ngspice_t24.txt'
file = fopen(filename, 'w');
fprintf (file, "Vs 1 0 0.0 ac 1.0 sin(0 1 1k)\nR1 1 2 %.11e\nR2 2 3 %.11e\nR3 2 5 %.11e\nR4 0 5 %.11e\nR5 5 6 %.11e\nR6 9 7 %.11e\nR7 7 8 %.11e\nVaux 0 9 0V\nHd 5 8 Vaux %.11e\nGb 6 3 (2,5) %.11e\nC1 6 8 %.11e\n\n.op\n\n.ic v(6) = %.11e v(8) = 0\n\n.end", R1, R2, R3, R4, R5, R6, R7, Kd, Kb, C, Vx);
fflush(filename);
fclose(filename);
