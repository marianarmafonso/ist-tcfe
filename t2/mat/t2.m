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
Ix = -(V6t0 - V5t0)*G(5) - Kb*(V2t0 - V5t0);
Req = abs(Vx/Ix);
tau = Req * C;
save("-ascii","../doc/tabela2.tex","V2t0","V3t0","V5t0","V6t0","V7t0","V8t0","Vx", "Ix", "Req", "tau");

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
Vf6 = V_amplitude(5) * sin(w*t) ; %solução forçada
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
xlabel('t[ms]');
ylabel('V');
title('Voltages as functions of time: v_s in grey and v_6 in black (Volts)');
print (hg, "../doc/v6.eps");

% ponto 6
numer= [0, 1];
denom= [Req*C,1];
numer2= [0,1];
denom2= [0, 1];
vc = tf(numer, denom); 
vs = tf(numer2,denom2);
w= logspace(-1, 6, 200); 
figure
bode(vc, vs, w);
print("lab2_bode.png", "-dpng");