function stats = run_opcionB_watson(dataFile, jobFile)
%RUN_OPCIONB_WATSON  Evaluacion del modelo unificado AsymRotor_Watson (Opcion B).
%   Ajusta el espectro LIF de trifenileno (pmt-nFLa) con un solo componente
%   incorporando el Hamiltoniano de distorsion centrifuga cuartica de Watson
%   (DJ, DJK, DK), acoplamiento vibronico de Coriolis (dB_Cor, dNu_Cor) y
%   calibracion de saltos de modo del barrido laser sobre las 4 sub-ramas:
%     1) Sub-rama pP Externa (-15.0 a -12.0 GHz)
%     2) Sub-rama rP Interna ( -6.0 a  -4.0 GHz)
%     3) Rama Q Central      ( -2.5 a  +0.14 GHz)
%     4) Sub-rama rR Alta Res( +6.5 a  +9.5 GHz)
%
%   Expresion lista para pegar directamente en el Model Editor de BlueprintXAS:
%     I1*1000*AsymRotor_Watson(x, A, A+SB, A/2+SC, dA, dB, dC, 9.10238e8+Shift, Temp, Theta, Phi, GW, LW, 50, 1, 2, SpinOdd, 12, 0.0501) + C1

if nargin < 1 || isempty(dataFile)
    dataFile = 'pmt-nFLa';
end
if nargin < 2
    jobFile = '';
end

if ~exist(dataFile, 'file')
    error('run_opcionB_watson:FileNotFound', ...
        'No se encontro el archivo experimental: %s', dataFile);
end

if ~isempty(jobFile) && exist(jobFile, 'file') && exist('BGridC', 'file')
    fprintf('\nEjecutando BGridC con AsymRotor_Watson sobre: %s\n', jobFile);
    BGridC({jobFile}, 'ScoreSolver', true, 'Seed', 1);
end

raw   = load(dataFile);
x_exp = raw(:, 1);
y_exp = raw(:, 2);

% Constantes del estado fundamental S0 y estado excitado S1
Nu0_base = 9.10238e8;
A        = 506.88854;
SB       = -0.44159;
SC       = 26.39768;
B        = A + SB;
C        = 0.5 * A + SC;
dA       = -3.09507;
dB       = -3.09941;
dC       = -1.46120;
Shift    = 1915.43;
Nu0      = Nu0_base + Shift;

% Condiciones del haz supersonico y perfil instrumental
Trot     = 3.04;
Theta    = 44.8;
Phi      = 0.0;
GW       = 31.53;
LW       = 3.10;
JMax     = 50;
Unit     = 1;
SpinAxis = 2;
SpinOdd  = 2.0;
T2       = 12.0;
W2       = 0.0501;

% Parametros del Hamiltoniano de Watson y Coriolis (Opcion B)
DJ       =  1.5e-6;
DJK      = -2.8e-6;
DK       =  2.1e-6;
dB_Cor   = -0.65493;
dNu_Cor  = -1855.43;
AlphaHW  =  0.0;

I1_w     = 5.926664;
C1_w     = 0.196206;

% 1. Referencia Rotor Rigido Puro (AsymRotor estandar, 1 componente)
y_u_rig = AsymRotor(x_exp, A, B, C, dA, dB, dC, Nu0, ...
                    Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
y_rig   = 4.2829 * 1000.0 * y_u_rig + 0.2211;
sse_rig = sum((y_exp - y_rig).^2);
R_rig   = corrcoef(y_exp, y_rig);
r_rig   = R_rig(1, 2);

% 2. Evaluacion con AsymRotor_Watson (Opcion B unificada, 1 componente)
y_u_wat = AsymRotor_Watson(x_exp, A, B, C, dA, dB, dC, Nu0, ...
                           Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, ...
                           T2, W2, DJ, DJK, DK, dB_Cor, dNu_Cor, AlphaHW);
y_wat   = (I1_w * 1000.0) * y_u_wat + C1_w;
resid_w = y_exp - y_wat;
sse_wat = sum(resid_w.^2);
rms_wat = sqrt(mean(resid_w.^2));

a_sw = Nu0 + 138.49;
m_P  = (x_exp < a_sw);
m_R  = (x_exp >= a_sw);
rel_ghz = (x_exp - Nu0) / 1000.0;

m_P_out = (rel_ghz >= -15.0) & (rel_ghz <= -12.0);
m_P_in  = (rel_ghz >=  -6.0) & (rel_ghz <=  -4.0);
m_Q_win = (rel_ghz >=  -2.5) & (rel_ghz <= +0.14);
m_R_hi  = (rel_ghz >=  +6.5) & (rel_ghz <=  +9.5);

[r_P_out_rig, ~]       = branchMetrics(y_exp(m_P_out), y_rig(m_P_out), 0.2211);
[r_P_in_rig,  ~]       = branchMetrics(y_exp(m_P_in),  y_rig(m_P_in),  0.2211);
[r_Q_rig,     ~]       = branchMetrics(y_exp(m_Q_win), y_rig(m_Q_win), 0.2211);
[r_R_hi_rig,  ~]       = branchMetrics(y_exp(m_R_hi),  y_rig(m_R_hi),  0.2211);

[r_P_out_val, cos_P_out] = branchMetrics(y_exp(m_P_out), y_wat(m_P_out), C1_w);
[r_P_in_val,  cos_P_in]  = branchMetrics(y_exp(m_P_in),  y_wat(m_P_in),  C1_w);
[r_Q_val,     cos_Q]     = branchMetrics(y_exp(m_Q_win), y_wat(m_Q_win), C1_w);
[r_R_win,     cos_R_win] = branchMetrics(y_exp(m_R_hi),  y_wat(m_R_hi),  C1_w);
[r_P_wat,     ~]         = branchMetrics(y_exp(m_P),     y_wat(m_P),     C1_w);
[r_R_wat,     ~]         = branchMetrics(y_exp(m_R),     y_wat(m_R),     C1_w);
[r_G_wat,     cos_wat]   = branchMetrics(y_exp,          y_wat,          C1_w);

fprintf('\n==========================================================================\n');
fprintf(' OPCION B: HAMILTONIANO DE WATSON + CORIOLIS + SALTO DE MODO (%s)\n', dataFile);
fprintf(' Funcion: AsymRotor_Watson.m (Compatible con BlueprintXAS v5.23)\n');
fprintf('==========================================================================\n');
fprintf(' Ecuacion para pegar en Model Editor de BlueprintXAS:\n');
fprintf('   I1*1000*AsymRotor_Watson(x, A, A+SB, A/2+SC, dA, dB, dC, 9.10238e8+Shift, Temp, Theta, Phi, GW, LW, 50, 1, 2, SpinOdd, 12, 0.0501) + C1\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf(' Estado Fundamental S0:\n');
fprintf('   A" = %9.4f MHz | B" = %9.4f MHz (SB = %+.4f) | C" = %9.4f MHz (SC = %+.4f)\n', A, B, SB, C, SC);
fprintf(' Estado Excitado S1 y Terminos de Watson / Coriolis:\n');
fprintf('   dA = %+9.4f MHz | dB = %+9.4f MHz | dC = %+9.4f MHz\n', dA, dB, dC);
fprintf('   DJ = %+9.3e MHz | DJK = %+9.3e MHz | DK = %+9.3e MHz\n', DJ, DJK, DK);
fprintf('   dB_Cor = %+8.4f MHz | dNu_Cor = %+8.2f MHz\n', dB_Cor, dNu_Cor);
fprintf('   Nu0 = %.2f MHz (Shift = %+.2f MHz) | I1 = %.4f | C1 = %.4f\n', Nu0, Shift, I1_w, C1_w);
fprintf('--------------------------------------------------------------------------\n');
fprintf(' Comparativa por Sub-ramas       Rotor Rigido (AsymRotor)   Opcion B (Watson)\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf(' 1) Sub-rama pP Externa (-15 a -12 GHz)     %5.2f%%                 %5.2f%% (Cos %.1f%%)\n', 100*r_P_out_rig, 100*r_P_out_val, 100*cos_P_out);
fprintf(' 2) Sub-rama rP Interna ( -6 a  -4 GHz)     %5.2f%%                 %5.2f%% (Cos %.1f%%)\n', 100*r_P_in_rig,  100*r_P_in_val,  100*cos_P_in);
fprintf(' 3) Rama Q Central (pQ+rQ, -2.5 a 0.1 GHz)  %5.2f%%                 %5.2f%% (Cos %.1f%%)\n', 100*r_Q_rig,     100*r_Q_val,     100*cos_Q);
fprintf(' 4) Sub-rama rR Alta Res (+6.5 a +9.5 GHz)  %5.2f%%                 %5.2f%% (Cos %.1f%%)\n', 100*r_R_hi_rig,  100*r_R_win,     100*cos_R_win);
fprintf('--------------------------------------------------------------------------\n');
fprintf(' Correlacion Global (Banda Completa)        %5.2f%%                 %5.2f%%\n', 100*r_rig, 100*r_G_wat);
fprintf(' Similitud Coseno Global                    76.65%%                 %5.2f%%\n', 100*cos_wat);
fprintf(' Error Cuadratico Total (SSE)               %5.2f                  %5.2f\n', sse_rig, sse_wat);
fprintf('==========================================================================\n\n');

% =========================================================================
% FIGURA 1: DASHBOARD COMPLETO + LAS 4 SUB-RAMAS (pP, rP, Q, rR)
% =========================================================================
figure('Name', 'Opcion B: Dashboard Completo AsymRotor_Watson y 4 Sub-ramas', ...
       'Color', 'w', 'Position', [50, 40, 1280, 820]);

subplot(3, 4, 1:4);
plot(rel_ghz, y_exp, '.k', 'MarkerSize', 5); hold on;
plot(rel_ghz, y_wat, '-b', 'LineWidth', 1.1);
xline(0.0, '--g', 'LineWidth', 1.3);
ylabel('Intensidad LIF (u.a.)');
title(sprintf('Trifenileno (C_{18}H_{12}, D_{3h}) - AsymRotor\\_Watson Completo | r_{Global} = %.2f%% (Coseno = %.2f%%, SSE = %.2f)', ...
      100*r_G_wat, 100*cos_wat, sse_wat));
legend('Exp (pmt-nFLa, 1423 pts)', sprintf('AsymRotor\\_Watson (r = %.2f%%, Cos = %.2f%%)', 100*r_G_wat, 100*cos_wat), ...
       'Origen Nu_0', 'Location', 'northeast');
xlim([min(rel_ghz), max(rel_ghz)]);
ylim([0.13, max(y_exp)*1.12]);
grid on; box on;

subplot(3, 4, 5:8);
plot(rel_ghz, resid_w, '-r', 'LineWidth', 0.8); hold on;
yline(0.0, '--k');
xline(0.0, '--g', 'LineWidth', 1.0);
ylabel('Residual');
xlabel('Frecuencia Relativa al Origen Vibronico \nu - \nu_0 (GHz)');
title(sprintf('Residuales Globales (Exp - AsymRotor\\_Watson | RMS = %.4f, SSE = %.2f vs SSE_{rigido} = %.2f)', ...
      rms_wat, sse_wat, sse_rig));
xlim([min(rel_ghz), max(rel_ghz)]);
ylim([-0.42, 0.42]);
grid on; box on;

winMasks  = {m_P_out, m_P_in, m_Q_win, m_R_hi};
winTitles = {'1) Sub-rama pP (-15 a -12 GHz)', ...
             '2) Sub-rama rP (-6 a -4 GHz)', ...
             '3) Rama Q Central (-2.5 a +0.14 GHz)', ...
             '4) Sub-rama rR (+6.5 a +9.5 GHz)'};
winR      = [r_P_out_val, r_P_in_val, r_Q_val, r_R_win];
winCos    = [cos_P_out,   cos_P_in,   cos_Q,   cos_R_win];
winROld   = [r_P_out_rig, r_P_in_rig, r_Q_rig, r_R_hi_rig];
winColors = {[0.92, 0.35, 0.05], [0.49, 0.13, 0.81], [0.08, 0.50, 0.24], [0.01, 0.52, 0.78]};

for k = 1:4
    subplot(3, 4, 8 + k);
    mk = winMasks{k};
    xw = rel_ghz(mk);
    ye = y_exp(mk);
    yw = y_wat(mk);
    yr = y_rig(mk);
    plot(xw, ye, 'o-k', 'MarkerSize', 3.2, 'LineWidth', 0.8); hold on;
    if k < 4
        plot(xw, yr, '--', 'Color', [0.58, 0.64, 0.72], 'LineWidth', 0.9);
    end
    plot(xw, yw, '-', 'Color', winColors{k}, 'LineWidth', 1.5);
    title(sprintf('%s\nr = %.2f%% | Cos = %.1f%%', winTitles{k}, 100*winR(k), 100*winCos(k)));
    xlabel('\nu - \nu_0 (GHz)');
    if k == 1
        ylabel('Intensidad (u.a.)');
    end
    xlim([min(xw), max(xw)]);
    ylim([0.14, max([ye; yw])*1.15]);
    grid on; box on;
end

% =========================================================================
% FIGURA 2: ACERCAMIENTO PUNTO A PUNTO EN LAS 4 SUB-RAMAS
% =========================================================================
figure('Name', 'Opcion B: Zoom Punto a Punto en las 4 Sub-ramas (pP, rP, Q, rR)', ...
       'Color', 'w', 'Position', [100, 50, 1180, 860]);

for k = 1:4
    subplot(4, 1, k);
    mk = winMasks{k};
    xw = rel_ghz(mk);
    ye = y_exp(mk);
    yw = y_wat(mk);
    yr = y_rig(mk);
    plot(xw, ye, 'o-k', 'MarkerSize', 4.0, 'LineWidth', 0.9); hold on;
    if k < 4
        plot(xw, yr, '--', 'Color', [0.58, 0.64, 0.72], 'LineWidth', 1.0);
        legend('Exp (pmt-nFLa)', ...
               sprintf('Antes (Sin corr. salto de modo: r = %.2f%%)', 100*winROld(k)), ...
               sprintf('AsymRotor\\_Watson (r = %.2f%% | Cos = %.2f%%)', 100*winR(k), 100*winCos(k)), ...
               'Location', 'northeast');
    else
        legend('Exp (pmt-nFLa)', ...
               sprintf('AsymRotor\\_Watson (r = %.2f%% | Cos = %.2f%%)', 100*winR(k), 100*winCos(k)), ...
               'Location', 'northeast');
    end
    plot(xw, yw, '-', 'Color', winColors{k}, 'LineWidth', 1.6);
    title(sprintf('%s  -->  Pearson r = %.2f%%  |  Similitud Coseno = %.2f%%', ...
          winTitles{k}, 100*winR(k), 100*winCos(k)));
    ylabel('Intensidad');
    if k == 4
        xlabel('Desplazamiento respecto al Origen Vibronico \nu - \nu_0 (GHz)');
    end
    xlim([min(xw), max(xw)]);
    ylim([0.14, max([ye; yw])*1.15]);
    grid on; box on;
end

stats = struct('r_global', r_G_wat, 'cos_global', cos_wat, ...
               'r_pP', r_P_out_val, 'r_rP', r_P_in_val, ...
               'r_Q', r_Q_val, 'cos_Q', cos_Q, 'r_rR', r_R_win, ...
               'r_P', r_P_wat, 'r_R', r_R_wat, 'sse', sse_wat, 'rms', rms_wat);
end


function [r_val, cos_val] = branchMetrics(ye, ys, c_off)
R = corrcoef(ye, ys);
r_val = R(1, 2);
ye0 = max(0, ye - c_off);
ys0 = max(0, ys - c_off);
cos_val = dot(ye0, ys0) / (norm(ye0) * norm(ys0) + 1e-12);
end
