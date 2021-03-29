R= [1003.32071212, 2044.60853047, 3082.91730437, 4160.61678649, 3040.22345043, 2067.11403452, 1033.02701196];
for i=1: 7
    G(i)=1/R(i);
    i=i+1;
end
Va= 5.13988034104; 
Id= 0.00102475824097;
Kb= 0.0070544535009;
Kc= 8161.13797582;
Z=0;
O=1;

% Método dos nós
A= [Z, -G(6), G(6)+G(7), Z, Z, Z, Z;
    Z, Z, Z, G(5), Z, Kb, -G(5)-Kb;
    Z, Z, Z, Z, G(2), -G(2)-Kb, Kb;
    -G(1), Z, Z, Z, -G(2), G(1)+G(2)+G(3), -G(3);
    O, -O, Z, Z, Z, Z, Z;
    Z, Kc*G(6), -Kc*G(6), Z, Z, Z, -O;
    G(1), G(4), G(7), Z, Z, -G(1), -G(4)]
B= [Z; Id; Z; Z; Va; Z; Z];
V=A\B;
V1 = V(1); 
V2 = V(2); 
V3 =V(3); 
V4 =V(4); 
V5 =V(5); 
V6 =V(6); 
V7 =V(7)

% Método das malhas
C= [R(4), Z, R(4)+R(6)+R(7)-Kc, Z;
    R(1)+R(3)+R(4), R(3), R(4), Z;
    Kb*R(3), -O+Kb*R(3), Z, Z;
    Z, Z, Z, O]
D=[Z; Va; Z; Id]
I=C\D;
I1 = I(1); 
I2 = I(2); 
I3 =I(3); 
I4 = I(4); 
% Criar tabela nós
save("-ascii","../doc/tabelaV.tex","V1", "V2", "V3", "V4", "V5", "V6", "V7");
save("-ascii","../doc/tabelaM.tex","I1", "I2", "I3", "I4");

Vc = Kc*I(3);
Ib = I(2);
IdM = I(4);
R1 = -I(1);
R2 = I(2);
R3 = I(1) + I(2);
R4 = -I(1) - I(3);
R5 = I(4) - I(2);
R6 = I(3);
R7 = I(3);
Vb = Ib/Kb;
V7M = Vc;
V6M = Vb + V7;
V5M = V6 + R(2)*I(2);
V4M = V7 + R(5)*(I(4)-I(2));
V3M = R(7)*I(3);
V2M = V3 + R(6)*I(3);
V1M = V2 + Va;

save("-ascii","../doc/tabelaM2.tex","Ib", "IdM", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "V1M", "V2M", "V3M", "V4M", "V5M", "V6M", "V7M");

IbN = Kb*(V(6) - V(7));
R1N = (V(6) - V(1))/R(1);
R2N = (V(5) - V(6))/R(2);
R3N = (V(6) - V(7))/R(3);
R4N = (V(2) - V(7))/R(4);
R5N = (V(4) - V(7))/R(5);
R6N = (V(2) - V(3))/R(6);
R7N = V(3)/R(7);

save("-ascii","../doc/tabelaV2.tex","IbN", "R1N", "R2N", "R3N", "R4N", "R5N", "R6N", "R7N");

