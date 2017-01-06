clear all; close all
%%
optimInfo.title = 'Função de referencia';
objFctHandle = @objetivo_allard_limp_funcao_dupla;
%cd('C:\Users\LVA 107\Documents\Dissertação\Capítulos\3 _ Caracterização inversa\Caracterização Tubo Grande\La de Rocha 8\Caracterização Inversa La de Rocha\GE')
cd('GE')

%%
paramDefCell = {'parameter1', [5000 300000], 10     % RESISITIVDADE AO FLUXO
                'parameter2', [0.90   0.99], 0.01   % POROSIDADE
                'parameter3', [1 4], 0.01           % TORTUOSIDADE
                'parameter4', [10e-6 500e-6], 1e-6  % COMRPIMENTO CARACTERÍSTICO VISCOSO
                'parameter5', [10e-6 500e-6], 1e-6  % COMPRIMENTO CARACTERÍSTICO TÉRMICO
                'parameter6', [64 64], 0.01 };      % DENSIDADE

objFctSettings = {};
objFctParams = [];
DEParams = getdefaultparams;
DEParams.saveHistory = 1;       % Se "0" nao salvar histórico se "1" salvar histórico
DEParams.NP = 100;              % TAMANHO DA POPULAÇÃO
DEParams.CR = 1;                % FATOR CROSSOVER - probabilidade de recombinação
DEParams.F = 0.935;             % FATOR MUTAÇÃO
% DEParams.feedSlaveProc  = 0;
DEParams.validChkHandle = @restricao;               % FUNÇÃO RESTRIÇÃO (Comprimento característico térmico igual ou maior que o comprimento caracteístico viscoso
DEParams.maxiter  = 150;        % NÚMERO MÁXIMO DE ITERAÇÕES
DEParams.maxtime  = inf;  % in seconds (tempo máximo não limitado "inf" = infito)
% DEParams.maxclock = [];
% DEParams.infoIterations = 0;
% DEParams.infoPeriod     = 10;  % in seconds
emailParams = [];
setrandomseed(1);

% start differential evolution
[bestmem, bestval, bestFctParams, nrOfIterations, resultFileName] = differentialevolution(...
  DEParams, paramDefCell, objFctHandle, objFctSettings, objFctParams, emailParams, optimInfo); 
disp(' ');
disp('Best parameter set returned by function differentialevolution:');
disp(bestFctParams);

% continue optimization by loading result file
% if DEParams.saveHistory
%   disp(' ');
%   disp(textwrap2(sprintf(...
%     'Now continuing optimization by loading result file %s.', resultFileName)));
%   disp(' ');
%   DEParams.maxiter  = 100;
%   DEParams.maxtime  = inf;  % in seconds
%   [bestmem, bestval, bestFctParams] = differentialevolution(...
%     DEParams, paramDefCell, objFctHandle, objFctSettings, objFctParams, emailParams, optimInfo, ...
%     resultFileName); 
%   disp(' ');
%   disp('Best parameter set returned by function differentialevolution:');
%   disp(bestFctParams);
% end

%%
%%
%%
% close all
Obj = bestval
cd('..')
load('A_exp')
load('Z_exp')
A_simples_exp = transpose(A_1);
A_duplo_exp = transpose(A_4);
f = 0:1:2048;

[Z_inv_simples A_inv_simples] = Allard_rigido(0.025,2048,bestmem(1),bestmem(2),bestmem(3),bestmem(4),bestmem(5));
[Z_inv_duplo A_inv_duplo] = Allard_rigido(2*0.025,2048,bestmem(1),bestmem(2),bestmem(3),bestmem(4),bestmem(5));

figure(20)
plot(f,A_simples_exp,'--k','linewidth',4);hold on; grid on
plot(f,A_duplo_exp,'k','linewidth',4);hold on; grid on
plot(f,A_inv_simples,'--b','linewidth',4); hold on; grid on
plot(f,A_inv_duplo,'b','linewidth',4); hold on; grid on
xlim([0 1800]); ylim([0 1])
legend('Experimento Simples','Experimento Duplo','Inverso Simples','Inverso Duplo')

%%
figure(21)
plot(f,real(Z_1)/413,'--k',f,imag(Z_1)/413,'--k','linewidth',3); hold on; grid on
plot(f,real(Z_4)/413,'k',f,imag(Z_4)/413,'k','linewidth',3); hold on; grid on
plot(f,real(Z_inv_simples)/413,'--b',f,imag(Z_inv_simples)/413,'--b','linewidth',3); hold on; grid on
plot(f,real(Z_inv_duplo)/413,'b',f,imag(Z_inv_duplo)/413,'b','linewidth',3); hold on; grid on

axis([ 100 1800 -4 4])