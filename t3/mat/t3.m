close all
clear all
pkg load symbolic;

% inicializar variáveis
C = 0.0004;
n_diodes = 18;
R_env= 400000;
R_reg= 57000;
V_primary = 230;
f = 50;
w = 2*pi*f;
Is = 1 * 10^(-14);
VT = 0.025; 
eta = 1;
n = 2.241016;
t = linspace(0, 0.2, 5000);

vc = zeros(1, length(t));
vo_enve = zeros(1, length(t));
vo_regac = zeros(1, length(t));

vo_reg = zeros(1, length(t));
vo_regdc = 0;
vo_regac = zeros(1, length(t));

% voltage regulator
ix = Is * (exp(12/(n_diodes*eta*VT))-1);%equacao do diodo
rd = eta * VT / (ix + Is); %resistencia incremental diodo
%V_r_reg = ix * R_reg; % voltagem nos terminais da resistência
%Vc = V_r_reg + 12; % lei das malhas -> variável devolve a voltagem que entra no circuito = voltagem (média) no condensador do envelope
%A_secondary = Vc;
A_secondary = V_primary/n;
for i=1:1:length(t)
    V_secondary (i) = A_secondary * cos(w*t(i));
end
%n = V_primary/A_secondary;

 save("-ascii","../doc/tabelaValues.tex","n_diodes", "n", "C", "R_env", "R_reg");
 
%envelope detector
    toff = (1/w) * atan(1/(w*R_env*C));
    
     for i= 1:1:length(t)
        if (V_secondary(i) > 0) %Diode on
            vc(i) = V_secondary(i);
        else %Diode off
            vc(i) = -V_secondary(i);
        end
     end

   
    % envelope
    for i=1:1:length(t)
        if t(i) < toff
            vo_enve(i) = vc(i);
        elseif A_secondary*abs(cos(w*toff))*exp(-(t(i)-toff)/(R_env*C)) > vc(i) % toff <= t < ton
            vo_enve(i) = A_secondary*abs(cos(w*toff))*exp(-(t(i)-toff)/(R_env*C)); %vo_exp(i)
        else % t >= ton
            toff = toff + 1/(2*f);
            %vo_exp (i) = A_secondary*cos(w*toff)*exp(-(t(i)-toff)/(R_env*C));
            vo_enve(i) = vc(i);
        end
    end
    
    average_voenve = mean(vo_enve);
    
    %Voltage regulator
    %Analysis dc
    if average_voenve >= 12
        vo_regdc = 12;
    else
        vo_regdc = average_voenve;
    end

    maximo_vo_reg = 0;
    minimo_vo_reg = 300000;
    
    %Incremental Analysis
    for i=1:1:length(t)
        if vo_enve (i) >= 12
            vo_regac(i) = ((n_diodes*rd)/(R_reg + n_diodes*rd))*(vo_enve(i) - average_voenve);
            vo_reg (i) = vo_regac (i)+ vo_regdc;
            if vo_reg(i) > maximo_vo_reg
                maximo_vo_reg = vo_reg(i);
            elseif vo_reg(i) < minimo_vo_reg
                minimo_vo_reg = vo_reg(i);
            end
        else
            vo_regac(i) = vo_enve(i) - average_voenve; 
            vo_reg (i)= vo_regac(i) + vo_regdc;
            if vo_reg(i) > maximo_vo_reg
                maximo_vo_reg = vo_reg(i);
            elseif vo_reg(i) < minimo_vo_reg
                minimo_vo_reg = vo_reg(i);
            end
        end
    end
   
    ripple = maximo_vo_reg - minimo_vo_reg ;
    average_voreg = mean(vo_reg);
    dev = mean(abs(average_voreg - 12));
    cost = R_env*10^-3 + R_reg*10^-3 + (4 + n_diodes)*0.1 + C*10^6;
    M = 1/(cost*(ripple + dev + 1*10^-6));
    
    save("-ascii","../doc/tabela1.tex","average_voreg", "ripple", "dev", "cost", "M");

hf =  figure();
plot (t*1000, vo_enve, "b");
legend("vo_{ENVELOPE}");
xlabel('t[ms]');
ylabel('V');
title('Voltage as function of time: vo_{ENVELOPE}');
print (hf, "../doc/plot1.eps");

hm =  figure();
plot (t*1000, vo_reg, "g");
legend("vo_{REGULATOR}");
xlabel('t[ms]');
ylabel('V');
title('Voltage as function of time: vo_{REGULATOR}');
print (hm, "../doc/plot2.eps");

hg =  figure();
plot (t*1000, vo_reg - 12, "r");
legend("Deviations");
xlabel('t[ms]');
ylabel('V');
title('Deviations from 12V');
print (hg, "../doc/plot3.eps");
