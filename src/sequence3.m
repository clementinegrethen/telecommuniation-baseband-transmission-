% séquence 3:

%% CHAINE DE REFERENCE SANS BRUIT:
% données
n=10000
close all;
Fe=24000 % en Hz
Rb=3000 % bits par seconde
Tb=1/Rb
Tb=1/Rb
nb_bit=5000
Nb=Fe/Rb %nb de point/bit
M=2
Ns=Nb*log2(M)
Ns1=Ns
n_oeil=1000
% mapping
message_binaire=randi([0,1],1,nb_bit);
Symboles = 2*message_binaire-1;
Suite_diracs = kron(Symboles, [1 zeros(1,Ns-1)]);

% Etude sans canal de propagation:
h = ones(1,Ns);
x = filter(h,1,Suite_diracs);
t = linspace(0,(Ns*nb_bit)/Rb,Ns*nb_bit);

% réception:
hr=ones(1,Ns);
xr=filter(hr,1,x);
%diagramme de l'oeil
oeil = reshape(xr, Ns, length(xr)/(Ns));
oeil = [0 oeil(end,1:(end-1)) ;  oeil(:,:)];
figure(14);
% réponse impulsionnelle globale
g = conv(h,hr);
figure(1);
plot(g);
% instant où n0 est optimal:
n0=Ns;


%échantillonnage du signal en sortie du filtre de reception (n0=Ns)
xe = xr(n0:Ns:end);
% taux d'erreur 
[bits_s] = (sign(xe)+1)/2;
check = find(message_binaire~= bits_s);
taux_erreur = length(check)/length(message_binaire)
%---------------------------------------------------------------------------

% On ajoute le bruit:

% réception
h = ones(1,Ns);
x = filter(h,1,Suite_diracs);
t = linspace(0,(Ns*nb_bit)/Rb,Ns*nb_bit);

% ajout du bruit
Px = mean(abs(x).^2); 
Eb_sur_N0 = 1;
sigma_carre = (Px*Ns)/(2*Eb_sur_N0);
x_avec_bruit = sqrt(sigma_carre)*randn(1,length(x)) +x;

% réception:
hr=ones(1,Ns);
xr=filter(hr,1,x_avec_bruit);

%diagramme de l'oeil:
%figure(3);%Première chaine à étudier:

plot(reshape(xr(Ns+1:end),Ns,length(xr(Ns+1:end))/Ns));

% taux erreur binaire: tracé et comparaison avec la théorie:

eps=5/100
max_db = 8
Px = mean(abs(x).^2);
teb = [];
EB = linspace(0,max_db, 9);
n0 = Ns;
i = 1;
for db = EB
    nb_err = 0;
    c = 0
    while nb_err < 1/(eps^2)
        c = c + 1;
        var_b = Px * Ns ./ (2*log2(M)*(10^(db/10)));
        bruit = sqrt(var_b)*randn(1,length(x));
        xb = x+bruit;
        zb = filter(hr,1,xb);
        zm = zb(n0:Ns:end);

        dm = 1*(zm > 0) + -1*(zm < 0);

        b1 = 1*(dm == 1) + 0*(dm == -1);

        erreur = message_binaire-b1;
        nb_err = nb_err + length(erreur( erreur ~= 0))
    end

    oeil = reshape(zb, Ns1, length(zb)/(Ns1));
    oeil = [0 oeil(end,end-n_oeil:(end-1)) ;  oeil(:,end-n_oeil:end)];

    subplot(3,3,i); plot(0:8, oeil(:,2:(end-1)));
    title(['Eb/N0 = ' , num2str(db) , 'dB'])
    i = i + 1 ;
    %Calcul de l'erreur
    taux_error = nb_err / (c * length(message_binaire));
    teb = [ teb taux_error]
end
figure(5);

semilogy(EB, teb);
hold on
grid on
teb_t = 1 -normcdf(sqrt(2*log2(M)*10.^(EB/10)))
teb_t_ref = teb
semilogy(EB, teb_t);
legend("TEB Simulé de la  chaine référence", "TEB Theorique")

%---------------------------------------------------------------------------

%% Première chaine à étudier:
% implantation de la chaine sans bruit
Fe=24000 % en Hz
Rb=3000 % bits par seconde
Tb=1/Rb
Tb=1/Rb
nb_bit=5000
Nb=Fe/Rb %nb de point/bit
M=2
Ns=Nb*log2(M)
% mapping
message_binaire=randi([0,1],1,nb_bit)
Symboles = 2*message_binaire-1
Suite_diracs = kron(Symboles, [1 zeros(1,Ns-1)])

% Etude sans canal de propagation:
h_1 = ones(1,Ns)
x_1 = filter(h,1,Suite_diracs)
t = linspace(0,(Ns*nb_bit)/Rb,Ns*nb_bit)

% réception:
hr= [ones(1,Ns/2) zeros(1,Ns/2)];

xr=filter(hr,1,x_1);
%diagramme de l'oeil:

oeil = reshape(xr, Ns, length(xr)/(Ns));
oeil = [0 oeil(end,1:(end-1)) ;  oeil(:,:)];
figure(16);
plot(0:8,oeil(:,2:end));
figure(9);
%échantillonnage du signal en sortie du filtre de reception (n0=Ns)
xe = xr(n0:Ns:end);
% taux d'erreur sans bruit:
[bits_s] = (sign(xe)+1)/2;
check = find(message_binaire~= bits_s);
taux_erreur_chaine_1_sans_bruit = length(check)/length(message_binaire)

%Ajout du bruit:
Px_1 = mean(abs(x_1).^2);
Eb_sur_N0 = 1;
sigma_carre_1 = (Px_1*Ns)/(2*Eb_sur_N0);
x_avec_bruit_1 = sqrt(sigma_carre)*randn(1,length(x_1)) +x_1;

% Filtrage de réception
hr_1(1:Ns) = ones(1,Ns);
hr_1(Ns/2+1:1:Ns) = 0;
xr_1 = filter(hr_1,1,x_avec_bruit_1);

% taux erreur binaire: tracé et comparaison avec la théorie:
max_db = 8;
Px = mean(abs(x_1).^2);
teb = [];
EB = linspace(0,max_db, 9);
n0 = Ns;
i=1
for db = EB
    nberr = 0;
    c = 0;
    while nberr < 1/(eps^2)
        c = c + 1;
        var_b = Px * Ns ./ (2*log2(M)*(10^(db/10)));
        bruit = sqrt(var_b)*randn(1,length(x_1));
        xb = x_1+bruit;
        zb = filter(hr,1,xb);
        zm = zb(n0:Ns:end);
        dm = 1*(zm > 0) + -1*(zm < 0);
        b1 = 1*(dm == 1) + 0*(dm == -1);
        erreur = message_binaire-b1;
        nberr = nberr + length(erreur( erreur ~= 0));
    end
    oeil = reshape(zb, Ns1, length(zb)/(Ns1));
    oeil = [0 oeil(end,end-n_oeil:(end-1)) ;  oeil(:,end-n_oeil:end)];

    subplot(3,3,i); plot(0:8, oeil(:,2:(end-1)));
    title(['Eb/N0 = ' , num2str(db) , 'dB'])
    i = i + 1 ;
   
    %Calcul de l'erreur
    taux_error = nberr / (c * length(message_binaire));
    teb = [ teb taux_error];
end
figure(21)
semilogy(EB, teb);
hold on
grid on
teb_t = 1 -normcdf(sqrt(log2(M)*10.^(EB/10)));
semilogy(EB, teb_t);
hold on
legend("TEB Simulé chaine 1", "TEB théorique")
figure(22)
semilogy(EB, teb);
hold on
grid on
semilogy(EB, teb_t_ref);
hold on
legend("TEB Simulé chaine 1", "TEB Référence")




%--------------------------------------------------------------------------
%% Deuxième chaine à étudier:
% implantation de la chaine sans bruit (première technique pour la partie sans bruit: elle m'aide à comprendre : pas forcement intéréssant pour le lecteur )
Fe=24000 % en Hz
Rb=3000 % bits par seconde
Tb=1/Rb
Tb=1/Rb
nb_bit=5000
Nb=Fe/Rb %nb de point/bit
M=4
Ns=Nb*log2(M);
% mapping
message_binaire=randi([0,1],1,nb_bit);
Symboles2 =(2*dec2bin(reshape(message_binaire,2,length(message_binaire)/2).')-3).';
Suite_diracs2 = kron(Symboles2, [1 zeros(1,Ns-1)]);

% Etude sans canal de propagation:
h_2 = ones(1,Ns);
x2 = filter(h,1,Suite_diracs2);
t = linspace(0,(Ns*nb_bit)/Rb,Ns*nb_bit);
% Filtrage de réception
hr2(1:Ns) = ones(1,Ns);
xr2 = filter(hr2,1,x2);




% n0 optimal = Ns

n0 = Ns;

% Echantillonnage

y2 = xr2(n0:Ns:end);
% symbolesDecides:

SymbolesDecides = zeros(length(y_2),1);
for i = 1:length(y_2)
    if y_2(i) > 2*Ns
        SymbolesDecides(i) = 3;
    elseif (y_2(i) >= 0) && (y_2(i) < 2*Ns)
        SymbolesDecides(i) = 1;
    elseif (y_2(i) < 0) && (y_2(i) > -2*Ns)
        SymbolesDecides(i) = -1;
    else
        SymbolesDecides(i) = -3;
    end
end



%% On ajoute le bruit:
% on ajoute le bruit (deuxième technique pour la partie sans bruit: elle est mieux )
Fe=24000 % en Hz
Te=1/Fe
Rb=3000 % bits par seconde
Tb=1/Rb
nb_bit=5000
Nb=Fe/Rb %nb de point/bit
M=4
Ns=Nb*log2(M);
%mapping
message_binaire=randi([0,1],1,nb_bit);
Symboles2 =(2*dec2bin(reshape(message_binaire,2,length(message_binaire)/2).')-3).';
Suite_diracs_2 = kron(Symboles2, [1 zeros(1,Ns-1)]);

t2 = linspace(0,Te*Ns*n, n*Ns/log2(M));
map = message_binaire*2-1;

map_2 = map(1:2:end).*(3-2*message_binaire(2:2:end));

s2 = kron(map_2,eye(1,Ns));

h2 = ones(1,Ns);
x2 = filter(h2,1,s2);

%Demodulation
hr = h2;
z2 = filter(hr,1,x2);
%diagramme de l'oeil:
oeil = reshape(z2, Ns, length(z2)/(Ns));
oeil = [0 oeil(end,1:(end-1)) ;  oeil(:,:)];
figure(11);
plot(0:16,oeil(:,2:(end-1)));
figure(12);
% n0 optimal
n0 = 16;
normz2 = z2(n0:Ns:end)*3/max(z2);


sign1 = sign(normz2);
sign2 = sign(normz2-2.*sign1);

sy2 = 3 * (sign1 == 1).* (sign2 == 1) + 1 * (sign1 == 1).* (sign2 == -1) + -1 * (sign1 == -1).* (sign2 == 1) + -3 * (sign1 == -1) .* (sign2 == -1);
demod1 = sign(kron(sy2, [1 0]));
demod2 = abs(kron(sy2, [0 1]));
b2 = (demod1 == 1) + (demod2  == 1);

erreur = message_binaire-b2;

%Demodulation avec bruit
Px = mean(abs(x2).^2);
teb_b = [];
teb_s = [];
Eb = linspace(0,max_db, 9);
n0 = Ns;
i = 1;
n_oeil=1000
for db = Eb
    nb_err_b = 0;
    nb_err_s = 0;
    compteur = 0;
    
    while nb_err_b < 1/(eps^2)
        ccompteur = ccompteur  + 1;
        var_b = Px * Ns ./ (2*log2(M)*(10^(db/10)));
        bruit = sqrt(var_b)*randn(1,length(x2));
        xb = x2+bruit;
        zb = filter(hr,1,xb);
        zm = zb(n0:Ns:end)*coeff;
        sign1 = sign(zm);
        sign2 = sign(zm-2.*sign1);
        sy2m = 3 * (sign1 == 1).* (sign2 == 1) + 1 * (sign1 == 1).* (sign2 == -1) + -1 * (sign1 == -1).* (sign2 == 1) + -3 * (sign1 == -1) .* (sign2 == -1);
        demod1 = sign(kron(sy2m, [1 0]));
        demod2 = abs(kron(sy2m, [0 1]));
        b2m = (demod1 == 1) + (demod2  == 1);
        erreur_b = message_binaire-b2m;
        erreur_s = sy2m - map2;
        nb_err_b = nb_err_b + sum(erreur_b ~= 0);
        nb_err_s = nb_err_s + sum(erreur_s ~= 0);
    end
    oeil = reshape(zb, Ns1, length(zb)/(Ns1));
    oeil = [0 oeil(end,end-n_oeil:(end-1)) ;  oeil(:,end-n_oeil:end)];

    subplot(3,3,i); plot(0:8, oeil(:,2:(end-1)));
    title(['Eb/N0 = ' , num2str(db) , 'dB'])
    i = i + 1 ;
    
   
   

    %Calcul de l'erreur
    taux_error_b = nb_err_b / (compteur * length(message_binaire));
    taux_error_s = nb_err_s / (compteur  * length(sy2m));
    teb_b = [ teb_b taux_error_b];
    teb_s = [ teb_s taux_error_s];
end
figure()

semilogy(Eb, teb_b);
hold on
grid on
Eb_No = log2(M)*10.^(Eb/10);
teb_t = (3/4)*(1-normcdf(sqrt((4/5)*power(10,([0:8]/10))),0,1))
semilogy(Eb, teb_t)
legend("TEB Simulé chaine 2", "TEB théorique ");
figure()
semilogy(Eb, teb_b);
hold on
grid on
semilogy(Eb, teb_t_ref)
legend("TEB Simulé chaine 2", "TEB théorique ");
figure()
teb_t_s=(1-normcdf(sqrt(4/5*Eb_No)))*3/2
semilogy(Eb,teb_s)
hold on 
grid on
semilogy(Eb,teb_t_s)
legend("TES simulé","TES de référence")

