% données de la séquence 2:
close all;
Fe=24000 % en Hz
Rb=3000 % bits par seconde
Tb=1/Rb
Tb=1/Rb
nb_bit=5000
Nb=Fe/Rb %nb de point/bit
M=2
Ns=Nb*log2(M);
% mapping
message_binaire=randi([0,1],1,nb_bit);
Symboles = 2*message_binaire-1;
Suite_diracs = kron(Symboles, [1 zeros(1,Ns-1)]);

%% Etude sans canal de propagation:
h = ones(1,Ns);
%modulateur 
x = filter(h,1,Suite_diracs);
t = linspace(0,(Ns*nb_bit)/Rb,Ns*nb_bit);
figure(1);
plot(t,x);
axis([0 (nb_bit-1)/Rb -1.5 1.5]);
xlabel('Temps (en s)');
ylabel('x(t)');
title('Signal transmis ');

% Calcul et tracé de la DSP
[DSP,f] = pwelch(x,[],[],[],Fe,'twosided');
figure(2);
semilogy(linspace(-Fe/2,Fe/2,length(DSP)),fftshift(DSP));
xlabel('Fréquence (en Hz)');
ylabel('DSP');
title('Densité spectrale de puissance du signal transmis ');

% filtre réception démodulateur :
hr=h;
xr=filter(hr,1,x);
figure(3);
plot(t,x);
axis([0 (nb_bit-1)/Rb -1.5 1.5]);
xlabel('Temps (en s)');
ylabel('x(t)');
title('Signal en sortie du filtre de réception ');

%échantillonnage du signal en sortie du filtre de reception (n0=Ns)
n0=Ns
xe = xr(n0:Ns:end);
figure(5); plot(xe);

% réponse impulsionnelle globale
g = conv(h,hr);

figure(22);
subplot(3,1,2), plot(g);
% instant où n0 est optimal:
n0=Ns;
%diagramme de l'oeil:
figure(4);
plot(reshape(xr(Ns+1:end),Ns,length(xr(Ns+1:end))/Ns));




% taux d'erreur 
[bits_s] = (sign(xe)+1)/2;
check = find(message_binaire~= bits_s);
taux_erreur = length(check)/length(message_binaire)

%échantillonnage du signal en sortie du filtre de reception (n0=3)
n0=3;
xe = xr(n0:Ns:end);
figure(6); plot(xe);

% taux d'erreur 
[bits_s] = (sign(xe)+1)/2;
check = find(message_binaire~= bits_s);
taux_erreur = length(check)/length(message_binaire)

%% Etude avec canal de propagation sans bruit:


Fe=24000 % en Hz
Rb=3000 % bits par seconde%échantillonnage du signal en sortie du filtre de reception (n0=Ns)
xe = xr(n0:Ns:end);
figure(5); plot(xe);

% taux d'erreur 
[bits_s] = (sign(xe)+1)/2;
check = find(message_binaire~= bits_s);
taux_erreur = length(check)/length(message_binaire)
Tb=1/Rb
Tb=1/Rb
nb_bit=5000
Nb=Fe/Rb %nb de point/bit
M=2
% réponse impulsionnelle globale
g = conv(h,hr);
figure(7);
plot(g);
Ns=Nb*log2(M);
message_binaire=randi([0,1],1,nb_bit);
Symboles = 2*message_binaire-1;
Suite_diracs = kron(Symboles, [1 zeros(1,Ns-1)]);

% pour une fréquence de coupure de 8000 Hz

fc=8000 %frequence de coupure en Hz
h = ones(1,Ns);
x = filter(h,1,Suite_diracs);

%filtre passe-bas pour le canal de propagantion
N=61
hc = (2*fc/Fe)*sinc(2*(fc/Fe)*[-(N-1)/2:(N-1)/2]);
x_filtre=filter(hc,1,[x zeros(1,(N-1)/2)]);
x_filtre = x_filtre((N-1)/2+1:end);

% réception:
hr=ones(1,Ns);
xr=filter(hr,1,x_filtre);

% réponse impulsionnelle globale% pour une fréquence de coupure de 8000 Hz

g = conv(h,hr);
g3=conv(g,hc);

figure(8);
plot(g3);

% tracé de |H(f)Hr(f)| et |Hc(f)|

figure(9);
plot(fftshift(abs(fft(g))));
hold on;
plot(fftshift(abs(fft(hc))));
title('Réponse en fréquence pour Fc=8000Hz');

%Diagramme de l'oeil:
figure(10);
plot(reshape(xr(Ns+1:end),Ns,length(xr(Ns+1:end))/Ns));

%échantillonnage du signal en sortie du filtre de reception (n0=Ns)
xe = xr(Ns:Ns:end);
figure(11); plot(xe);

% taux d'erreur 
[bits_s] = (sign(xe)+1)/2;
check = find(message_binaire~= bits_s);
taux_erreur8000 = length(check)/length(message_binaire)

% pour une fréquence de coupure de 1000 Hz:


fc=1000 %frequence de coupure en Hz
h = ones(1,Ns);
x = filter(h,1,Suite_diracs);

%filtre passe-bas pour le canal de propagantion
N=61
hc2 = (2*fc/Fe)*sinc(2*(fc/Fe)*[-(N-1)/2:(N-1)/2]);
x_filtre=filter(hc2,1,[x zeros(1,(N-1)/2)]);
x_filtre = x_filtre((N-1)/2+1:end);

% réception:
hr=ones(1,Ns);
xr=filter(hr,1,x_filtre);

% réponse impulsionnelle globale% pour une% pour une fréquence de coupure de 8000 Hz fréquence de coupure de 1000 Hz

g1 = conv(h,hr);
g2=conv(g1,hc2);
figure(12);
plot(g2);

% tracé de |H(f)Hr(f)| et |Hc(f)|

figure(13);
plot(fftshift(abs(fft(g1))));
hold on;
plot(fftshift(abs(fft(hc))));
title('Réponse en fréquence pour Fc=1000Hz');


%Diagramme de l'oeil:
figure(14);
plot(reshape(xr(Ns+1:end),Ns,length(xr(Ns+1:end))/Ns));

%échantillonnage du signal en sortie du filtre de reception (n0=Ns)
xe = xr(Ns:Ns:end);
figure(15); plot(xe);

% taux d'erreur 
[bits_s] = (sign(xe)+1)/2;
check = find(message_binaire~= bits_s);
taux_erreur1000 = length(check)/length(message_binaire)






















