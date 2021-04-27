
%voltage regulator
ita = 2;
Vt = 0.025;
vd = 0.7;
Is = 1*(10^-14);
rd = (ita*Vt)/(Is*exp(vd/(ita*Vt)));
R = 15*(10^3);

