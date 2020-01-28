
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%     Programa para a solução do problema de Fluxo de Potência via       %
%                via método de Newton-Raphson                            %
%          (formulação convencional em coordenadas cartesianas)          %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



disp(' ')
disp('     ***************************************************')
disp('     *   Programa de Fluxo de Potência via  método     *')
disp('     *   de Newton-Raphson - coordenadas retangulares  *')
disp('     *           (formulação convencional)             *')
disp('     ***************************************************')
disp('     LABSPOT/UFSC                          Fev de 2008')
disp('     R. S. Salgado')
disp(' ')

caso = input(' nome do arquivo de dados (sem extensão):   ', 's');

if ~exist(caso)
  disp(' ')
  disp('   ********************************')
  disp('   *     ARQUIVO INEXISTENTE!     *')
  disp('   ********************************')
  disp(' ')
  return;
end

   format short e;
   format compact;

disp(' ')
disp('   Iteração   Max_Desbalanço       miv          ||dx||')
disp(' ========================================================')


   eval(caso);


%-----  especificação das constantes  -----

   Sbase = 100;             % potência aparente base
   radeg = 180/pi;
   tol_res = 1.0e-04;
   tol_mi = 1.0e-04;
   maxit = 50;
   j = sqrt(-1);
   fc = 1.0;                % fator de carga


   %-----  LEITURA DOS DADOS DO SISTEMA  -----
   
LT(7,:)=[];

nlt = size(LT, 1);
nb = size(barras, 1);

nsb=LT(:,2);
neb=LT(:,3);

no = LT(:, 1);
nd = LT(:, 2);
circ = LT(:,3);
r = LT(:, 4) / 100.0;
x = LT(:, 5) / 100.0;
b = LT(:, 6) / 200.0;
a = LT(:, 7);
amin = LT(:,8);
amax = LT(:,9);


% faca a renumeracao das barras e dos circuitos

AuxBarr = barras(:, 1);

for ij = 1:nb
   barras(ij, 1) = ij;
end 

for ik = 1 : nlt
    for jk = 1:nb
        if (no(ik) == AuxBarr(jk))
            nsb(ik) = jk;
        end
       if (nd(ik) == AuxBarr(jk))
           neb(ik) = jk; 
       end
   end 
end


barra = barras(:, 1);
tipo = barras(:, 2);
V = barras(:, 3);
Pg = barras(:, 5) / Sbase;
Qg = barras(:, 6) / Sbase;
Qgmin = barras(:, 7) / Sbase;
Qgmax = barras(:, 8) / Sbase;
Pd = fc * barras(:, 9) / Sbase;
Qd = fc * barras(:, 10) / Sbase;
shunt = barras(:, 11) / Sbase;
iAr  = barras(:, 12);

%-----  LEITURA DOS LIMITES de tensão e potência ativa gerada e custo de geração -----

Vmin = Lim(:, 2);
Vmax = Lim(:, 3);
Pgmin = Lim(:, 4) / Sbase;
Pgmax = Lim(:, 5) / Sbase;
bg = Lim(:, 6);
cg = Lim(:, 7);

for ibb=1:nb
     Lim(ij,1) = ibb;
end 
   
%-----  especificação dos conjuntos de barras  -----

   bpv = find(tipo == 1 | tipo == 4);
   npv = length(bpv);
   bpq = find(tipo == 0 | tipo == 3);
   npq = length(bpq);
   bf  = find(tipo == 2);
   bsf = find(tipo ~= 2);
   nsf = nb - 1;
   bge = find(tipo == 2 | tipo == 1 | tipo == 4);
   nge = length(bge);

%-----  formação da lista de adjacência e da matriz Y_barra -----

%-----  formação da lista de adjacência  -----

l = 1;
for i = 1 : nb
  no_temp = 0;
  nd_temp = 0;
  adj(i) = l;
  for k = 1 : nlt
    if (nsb(k) ~= no_temp) | (neb(k) ~= nd_temp)
      if nsb(k) == i
        b_adj(l) = neb(k);
        l = l + 1;
      end
      if neb(k) == i
        b_adj(l) = nsb(k);
        l = l + 1;
      end
    end
    no_temp = nsb(k);
    nd_temp = neb(k);
  end
end
adj(nb + 1) = l;


%-----  determinação da matriz Y_barra  -----%

B1_barra = zeros(nb, nb);


y = 1.0 ./x;

for l = 1 : nlt
    
  i = nsb(l);
  k = neb(l);

    B1_barra(i, i) = B1_barra(i, i) +  y(l);
    B1_barra(k, k) = B1_barra(k, k) + y(l);
    B1_barra(i, k) = B1_barra(i, k) - y(l);
    B1_barra(k, i) = B1_barra(i, k);
  end

  
B1_barra(1,:)=[];
B1_barra(:,1)=[];
B1_barra=B1_barra';

   % -----  cálculo das injeções de potência ----- %
   
P=Pg-Pd;
Pref=P(1);
P(1)=[];   

   % ------ atualização B_barra ---------- %
   
B1 = inv(B1_barra);

%% -----  solução do sistema linear  -----

delta=B1*P;

 % Conversão outputs
 
delta=delta*radeg;
delta=[0;delta];
P=[Pref;P];
Ptot = sum(P);


disp(' ')
fprintf( '\n==============================================================================');
fprintf( '\n|              Resultados das barras Fluxo Linearizado Contingênciado        |');
fprintf( '\n==============================================================================');
fprintf( '\n  Barra Tipo  abertura angular               Potência                           ');
fprintf('\n #                  graus                       P(MW)                                       ');
fprintf('\n------------------------------------------------------------------------------');
for i = 1:nb
    V(i)=1;
   fprintf('\n  %3.0f %3.0f          %6.2f                   %8.2f  \n',AuxBarr(i),  tipo(i),    delta(i),     P(i)');
end
fprintf('\n                                                             ');
fprintf( '\nPotência Total:                             %8.2f     ', Ptot);
fprintf('\n');
