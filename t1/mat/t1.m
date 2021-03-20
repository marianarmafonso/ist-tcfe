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
TN=[1, V(1); 2, V(2); 3, V(3); 4, V(4); 5, V(5); 6, V(6); 7, V(7)]

% Método das malhas
C= [R(4), Z, R(4)+R(6)+R(7)-Kc, Z;
    R(1)+R(3)+R(4), R(3), R(4), Z;
    Kb*R(3), -O+Kb*R(3), Z, Z;
    Z, Z, Z, O]
D=[Z; Va; Z; Id]
I=C\D;
TM=[1, I(1); 2, I(2); 3, I(3); 4, I(4)]

% Criar tabela nós
T = table(TN, 'VariableNames', {'Node Number', 'Node Voltage'});
L = table(TM, 'VariableNames', {'Loop Number', 'Loop Current'});
