%--------------------------------------------------------------------------
% Detta skript simulerar den cirkadiska rytmen med en stokastisk 
% metod, nämligen Gillespies algoritm. Vi tittar på samverkan     
% mellan proteinerna aktivator A och repressor R och plottar      
% detta. Eftersom det blir ett stort antal tidssteg sparar vi    
% endast data i % vart 100:e steg för plottning.      
%
%--------------------------------------------------------------------------

clear all;

% Totala tiden för simuleringen 
T= 200;

% Parametervärden
alfaA = 50; 
alfapA = 500; 
alfaR = 0.01; 
alfapR = 50; 
betaA = 50; 
betaR = 5; 
tetaA = 50; 
tetaR = 100;
gammaA = 1; 
gammaR = 1; 
gammaC = 2;
deltaMR = 0.5; 
deltaMA = 10; 
deltaA = 1; 
deltaR = input('Ge ett värde på den spontana nedbrytningen av R (deltaR): ');
 
% Spara alla parametrar i en vektor p 
p = [alfaA; alfapA; alfaR; alfapR; betaA; betaR; tetaA; tetaR; gammaA; gammaR; gammaC; deltaMR; deltaMA; deltaA; deltaR]; 

% Skapa initiala tillståndsvektorn x.
A0 = 0; 
C0 = 0;
DA0 = 1; 
DAp0 = 0; 
DR0 = 1; 
DRp0 = 0; 
MA0 = 0; 
MR0 = 0; 
R0 = 0; 
 

x = [A0 C0 DA0 DAp0 DR0 DRp0 MA0 MR0 R0]; 

% Initiala tiden t0
t=0;

% Läs in matris innehållandes alla stökiometri-vektorer för alla 
% reaktioner r (har 18 st). Dessa finns lagrade i nr_vilar.m
nr_matris = nr_vilar(); 

% Vi vill spara tillståndsvektorer x från några tidssteg för plottning.
% Börja med att fylla pa initiala tillståndsvektorn x:  
x_matris(1,:)=x; 

% Vi vill även spara motsvarande tidpunkter.Dessa sparas i t_vektor.  
t_vektor(1)=t; 

% Skapa index-variabel (raknare) far att halla koll pa vilket        % tidssteg vi befinner oss i. 
i=1;  

% Skapa separat index-variabel (räknare) för indexering av de        % vektorer/matriser vi vill spara for plottning (x_matris och        % t_vektor). 
i_spara=2;  %Har redan sparat värden pa första platsen => i_spara=2

% TIDS-LOOP
while t<T   
    % Läs in vektor (prop_vilar()) innehållandes propensiteter. 
    w = prop_vilar(x,p); 
    
    % Summera alla propensiteter.
    a0= sum(w); 
       

    %------------
    % Slumpa fram tidssteg tau fran exponentialfördelning
    % m.h.a. INVERSE TRANFORM SAMPLING:
    %------------
    % a.) Skapa ett slumpmässigt tal från likformig fördelning           % U(0,1)      
    u=rand(1,1);
    % b.) Beräkna slumpmässigt tidssteg tau fran                                %  exponentialfördelning 
    %     genom: F(tau)=1-e^(-a0*tau)   =>  F_inv(u)=-log(1-u)/a0
    tau=-log(1-u)/a0; 
    
    
    %------------
    % Slumpa fram nästa reaktion m.h.a. INVERSE TRANFORM SAMPLING:
    %------------
    % a.) Använd samma slumpat u som ovan.
    % b.) Vilken reaktion r som skall ske är en diskret stokastisk
    %     variabel.
    %     Berakna kumulativa sannolikheter F(r) för alla reaktioner, 
    %     spara dessa i vektor F:
    F=cumsum(w)./a0; 
    % Hitta sedan den reaktion r som gör att kumulativa                 % sannolikheten F(r)>u
    r=find(F>u,1);      
    
    % Uppdatera tillståndsvektorn x med hjälp av stökiometri-vektorn for reaktion r. 
    x=x+ nr_matris(r,:); 
    
    % Uppdatera tidssteget
    t=t+tau;
    
    % Spara lösningar från vissa tidssteg för plottning:
    if (mod(i,100)==0)  % => var hundrade tidssteg
        t_vektor(i_spara)= t;   
        x_matris(i_spara,:)= x;   
        i_spara = i_spara+1;       
                              
    end
    
    % öka index (räknaren) för vilket tidssteg vi befinner oss i:
    i= i+1;   
end

% Plotta A och R mot de sparade tidsstegen:
figure; 
plot(t_vektor, x_matris(:,1));
hold on; 
plot(t_vektor,x_matris(:,9)); 
title('Stokastisk simulering av A och R i den cirkadiska rytmen')
xlabel('Tid (h)');
ylabel('Antal molekyler');
legend('A','R');
