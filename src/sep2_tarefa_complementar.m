
disp(' ')
disp('     ***************************************************')
disp('     *   Programa de Fluxo de Potência Linearizado via *')
disp('     *  Análise de Sensibilidade/Injeções Compensadoras   *')
disp('     *                      *')
disp('     ***************************************************')
disp('             Sistemas Elétricos de Potência II ')
disp('    ')
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
   
nlt = size(LT, 1);
nb = size(barras, 1);


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

lt_tr = find(amin < amax);
ntraf = length(lt_tr);

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

B_barra = zeros(nb, nb);


y = 1.0 ./ x;

for l = 1 : nlt
  i = nsb(l);
  k = neb(l);

    B_barra(i, i) = B_barra(i, i) + y(l);
    B_barra(k, k) = B_barra(k, k) + y(l);
    B_barra(i, k) = B_barra(i, k) - y(l);
    B_barra(k, i) = B_barra(i, k);

end
B_barra(1,:)=[];
B_barra(:,1)=[];
B_barra=B_barra';

  % -----  cálculo das injeções de potência ----- %
   
P=Pg-Pd;
Pref=P(1); 
P(1)=[];
   % atualização B_barra
   
B = inv(B_barra);

%-----  solução do sistema linear  -----%

  ang=B*P;
  P=B_barra*ang; 
   
%-----  cálculos complementares  -----

Ptot = sum(P);

tipo(bge) = 1;
tipo(bf) = 2;
bpv = find(tipo == 1);
npv = length(bpv);
bpq = find(tipo == 0 | tipo == 3);
npq = length(bpq);


% ----- impressão dos resultados ------- %

V=1;
I=P/V;
ang=ang*radeg;

%% Análise de Contingência %%


%% Calculo do incremento angular por analise de sensibilidade

% ekl=vetor com elementos nulos e outros elementos ek=1 e el=-1

 ekl=[0;0;-1;1;zeros(9,1)];
 
% Matriz W

Z0=inv(B_barra);

W=Z0*ekl;

%Delta Y
deltay=-1/(x(7));

%% --- Sensibilidade angulo --- %%

delta_ang=-(deltay*(ang(4)-ang(3)))/(1+(ekl')*W*deltay);
delta_ang=delta_ang*W;


%% ----- Injeções Compensadoras -------- %%

tkl=-deltay*(ang(4)-ang(3));
Zeq=Z0(3,3)+Z0(4,4)-2*Z0(4,3);
delta_P=(tkl)/(1-(Zeq*(1/x(7))));
delta_P=delta_P*ekl;
deltaang=Z0*(delta_P);



% Recolocação das variáveis da barra de folga
P=[Pref;P];
ang=[0;ang];
delta_ang=[0;delta_ang];
deltaang=[0;deltaang];
delta_P=[0;delta_P];

% Calculo dos ângulos pela sensibilidade e pela injeção compensadora
sol=ang+delta_ang;
sol2=ang+deltaang;
Ptot = sum(P);
delta_P=delta_P/radeg;

disp(' ')
fprintf( '\n==============================================================================');
fprintf( '\n|            Resultado Análise de Sensibilidade/Injeções Compensadoras       |');
fprintf( '\n==============================================================================');
fprintf( '\n  Barra Tipo  Ângulos Sensibilidade Ângulo Novo  Injeção  Ângulo Novo  Potência ');
fprintf('\n              Caso Base                                                         ');
fprintf('\n #              graus      graus        graus                graus       P(MW)  ');
fprintf('\n-------------------------------------------------------------------------------');
for i = 1:nb
    V(i)=1;
   fprintf('\n  %3.0f %3.0f     %6.2f     %6.2f       %6.2f     %6.2f     %6.2f   %8.2f  \n',AuxBarr(i),  tipo(i),    ang(i),    delta_ang(i),    sol(i),  delta_P(i),  sol2(i),  P(i)');
end
fprintf('\n                                                                                ');
fprintf( '\nPotência Total:                                                      %8.2f     ', Ptot);
fprintf('\n Caso Base');

