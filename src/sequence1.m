
%modulateur1
%Variables
Fe=24000 %fréquence d'échantillonnage en Hz
Rb=3000 %débit binaire en bits par seconde
Tb=1/Rb% durée binaire
nb_bit=5000
Nb=Fe/Rb %nb de point/bit
M=2% mapping binaire
Ns=Nb*log2(M);%
message_binaire=randi([0,1],1,nb_bit);
nb_bits=1000
% mapping
Symboles = 2*message_binaire-1;
Suite_diracs = kron(Symboles, [1 zeros(1,Ns-1)]);
% modulation
h = ones(1,Ns);
x = filter(h,1,Suite_diracs);

figure; plot(x);
axis([0 nb_bits-1 -1.5 1.5]);

% Calcul et tracé de la DSP
[DSP,f] = pwelch(x,[],[],[],Fe,'twosided');
figure;
semilogy(linspace(-Fe/2,Fe/2,length(DSP)),fftshift(DSP));
xlabel('Fréquence (en Hz)');
ylabel('DSP');
title('Densité spectrale de puissance du signal transmis (modulateur 1)');

%DSP théorique 
f=linspace(-Fe/2,Fe/2,length(DSP));
Rs=Rb
Ts=Tb
DSP_theorique=Ts*power(sinc(f*Ts),2)
figure;
semilogy(f,DSP_theorique);
hold on;
semilogy(linspace(-Fe/2,Fe/2,length(DSP)),fftshift(DSP));
xlabel('Fréquence (en Hz)');
title('Comparaison avec la DSP théorique (modulateur 1)');


%Modulateur 2:
M2=4
Ns2=Nb*log2(M2);
%mapping
message_binaire2=randi([0,1],nb_bit/2,2);
symboles2=2*bi2de(message_binaire2)-3;
suite_diracs2=kron(symboles2',[1 zeros(1,Ns2-1)]);
h2=ones(1,Ns2)
%modulateur
x=filter(h,1,suite_diracs2)
t = linspace(0,(Ns2*nb_bit)/Rb,Ns2*nb_bit/2);
figure;
plot(t,x);
axis([0 (nb_bit-1)/Rb -1.5 1.5]);
xlabel('Temps (en s)');
ylabel('x(t)');
title('Signal transmis du modulateur 2');


% Calcul et tracé de la DSP
[DSP_Mod2,f] = pwelch(x,[],[],[],Fe,'twosided');
figure;
semilogy(linspace(-Fe/2,Fe/2,length(DSP_Mod2)),fftshift(DSP_Mod2));
xlabel('Fréquence (en Hz)');
ylabel('DSP');
title('Densité spectrale de puissance du signal transmis (modulateur 2)');

% DSP théorique
Rs=Rb/2
Ts=1/Rs
f=linspace(-Fe/2,Fe/2,length(DSP_Mod2));
DSP_theorique_MOD2=Ts*power(sinc(f*Ts/2),2);
figure;
semilogy(f,DSP_theorique_MOD2);
hold on;
semilogy(linspace(-Fe/2,Fe/2,length(DSP_Mod2)),fftshift(DSP_Mod2));
xlabel('Fréquence (en Hz)');
title('Comparaison avec la DSP théorique (modulateur 2)');



% Implantation du modulateur 3
Fe = 24000;
nb_bits = 1000;
%mapping
bits = randi([0,1],1,nb_bits);
Symboles = 2*bits-1;
Suite_diracs = kron(Symboles, [1 zeros(1,Ns-1)]);
h = rcosdesign(0.6,8,8);
%modulateur 
x = filter(h,1,Suite_diracs);

t = linspace(0,(Ns*nb_bits/2)/Rb,Ns*nb_bits);
figure;
plot(t,x);
axis([0 (nb_bits-1)/Rb -1.5 1.5]);
xlabel('Temps (en s)');
ylabel('x(t)');
title('Signal transmis du modulateur 3');



% Calcul et tracé de la DSP
[DSP_Mod3,f] = pwelch(x,[],[],[],Fe,'twosided');
figure;
semilogy(linspace(-Fe/2,Fe/2,length(DSP_Mod3)),fftshift(DSP_Mod3));
xlabel('Fréquence (en Hz)');
ylabel('DSP');
title('Densité spectrale de puissance du signal transmis (modulateur 3)');

% DSP théorique

f=linspace(-Fe/2,Fe/2,length(DSP_Mod3));
Rs=Rb
Ts=Tb
alpha=0.6;
var_3=var(x);
DSP_theorique_mod3=f;
v1=(1-alpha)/(2*Ts);
v2=(1+alpha)/(2*Ts);
for i=1:length(f)
    if  abs(f(i))<v1
        DSP_theorique_mod3(i)=var_3;
    elseif abs(f(i))<v2
        DSP_theorique_mod3(i)=Ts/2*(1+cos((pi*Ts/alpha)*(abs(f(i))-v1)))*(var_3/Ts);
    else 
        DSP_theorique_mod3(i)=0;
    end
end

figure;
semilogy(f,DSP_theorique_mod3/max(DSP_theorique_mod3));
hold on;
semilogy(linspace(-Fe/2,Fe/2,length(DSP_Mod3)),fftshift(DSP_Mod3/max(DSP_Mod3)));
xlabel('Fréquence (en Hz)');
title('Comparaison avec la DSP théorique (modulateur 3)');







%% Tracé des 3 DSP sur une figure:
figure
semilogy(linspace(-Fe/2,Fe/2,length(DSP)),fftshift(DSP));
hold on;

semilogy(linspace(-Fe/2,Fe/2,length(DSP_Mod2)),fftshift(DSP_Mod2));
hold on;
semilogy(linspace(-Fe/2,Fe/2,length(DSP_Mod3)),fftshift(DSP_Mod3));
hold on;

xlabel('Fréquence (en Hz)');
legend("modulateur1", "modulateur2","modulateur3")
title('tracé des 3 DSP');




