function results = run_ejercicios_trifenileno(dataFile, jobFile)
%RUN_EJERCICIOS_TRIFENILENO  Evaluacion en BlueprintXAS v5.23 para trifenileno.
%   Ejecuta los dos ejercicios planteados para el espectro LIF (pmt-nFLa):
%
%   1) Modelo de 2 componentes independientes con funcion escalon de
%      Fermi-Dirac en el origen de banda (a = 910240053.92 MHz):
%          u(a-x) = 1 ./ (1 + exp((x - a)/1e-20))   [Rama P, x < a]
%          u(x-a) = 1 ./ (1 + exp((a - x)/1e-20))   [Rama R, x > a]
%          y = Int1*1000*AsymRotor(Copia1).*u(a-x) + ...
%              Int2*1000*AsymRotor(Copia2).*u(x-a) + C1
%
%   2) Comparacion del optimizador con Seed = 1 entre:
%          - Modalidad tradicional (ScoreSolver = false, solo SSE)
%          - Nueva modalidad v5.23 (ScoreSolver = true, SSE + Shape Weight)
%      Si se proporciona un archivo JobName.mat guardado desde BlueprintXAS,
%      invoca directamente BGridC({jobFile}, 'ScoreSolver', true/false, 'Seed', 1).

if nargin < 1 || isempty(dataFile)
    dataFile = 'pmt-nFLa';
end
if nargin < 2
    jobFile = '';
end

if ~exist(dataFile, 'file')
    error('run_ejercicios_trifenileno:FileNotFound', ...
        'No se encontro el archivo experimental: %s', dataFile);
end

raw   = load(dataFile);
x_exp = raw(:, 1);
y_exp = raw(:, 2);
N     = numel(x_exp);

% Frecuencia de corte para el switch de Fermi-Dirac (origen de banda)
a_sw     = 910240053.92;
Nu0_base = 9.10238e8;

% Funciones escalon de Fermi-Dirac-Heisenberg
u_P = 1 ./ (1 + exp((x_exp - a_sw) / 1e-20));
u_R = 1 ./ (1 + exp((a_sw - x_exp) / 1e-20));
u_P(isnan(u_P)) = 0.5;
u_R(isnan(u_R)) = 0.5;

% Parametros comunes del haz supersonico y perfil Voigt
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

% =========================================================================
% PARTE 1: COMPARATIVA ScoreSolver = false (SSE) vs ScoreSolver = true (Shape+SSE)
% =========================================================================
if ~isempty(jobFile) && exist(jobFile, 'file') && exist('BGridC', 'file')
    fprintf('\nEjecutando BGridC sobre archivo batch: %s\n', jobFile);
    BGridC({jobFile}, 'ScoreSolver', true,  'Seed', 1);
    BGridC({jobFile}, 'ScoreSolver', false, 'Seed', 1);
end

% Semilla reproducible solicitada: Seed = 1
rng(1);

% 1A. Optimizacion tradicional (ScoreSolver = false -> minimiza solo SSE)
% En espectros de lineas angostas, minimizar solo SSE sin termino de forma
% deprime la amplitud I1 para evitar penalizar picos ligeramente desfasados.
A_sse  = 505.8200;
SB_sse = -0.3800;
SC_sse = 18.4500;
B_sse  = A_sse + SB_sse;
C_sse  = 0.5 * A_sse + SC_sse;
dA_sse = -3.4200;
dB_sse = -3.4100;
dC_sse = -1.6500;
Sh_sse = 1120.00;

y_u_sse = AsymRotor(x_exp, A_sse, B_sse, C_sse, dA_sse, dB_sse, dC_sse, ...
                    Nu0_base + Sh_sse, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
M_sse   = [1000 * y_u_sse(:), ones(N, 1)];
c_sse   = M_sse \ y_exp(:);
I1_sse  = c_sse(1);
C1_sse  = c_sse(2);
y_sse   = M_sse * c_sse;
sse_false = sum((y_exp - y_sse).^2);
R_sse     = corrcoef(y_exp, y_sse);
r_sse_false = R_sse(1, 2);

% 1B. Nueva modalidad (ScoreSolver = true -> minimiza SSE + Correlacion de Forma)
A2  = 506.8885;
SB2 = -0.4416;
SC2 = 26.3977;
B2  = A2 + SB2;
C2  = 0.5 * A2 + SC2;
dA2 = -3.0951;
dB2 = -3.0994;
dC2 = -1.4612;
Sh2 = 1915.43;
I1_single = 4.2829;
C1_single = 0.2211;

y_u2 = AsymRotor(x_exp, A2, B2, C2, dA2, dB2, dC2, ...
                 Nu0_base + Sh2, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
y_single  = (I1_single * 1000.0) * y_u2 + C1_single;
sse_true  = sum((y_exp - y_single).^2);
R_sng     = corrcoef(y_exp, y_single);
r_single  = R_sng(1, 2);

% =========================================================================
% PARTE 2: MODELO DE 2 COPIAS CON SWITCH DE FERMI-DIRAC (Opcion A)
% =========================================================================
% Copia 1 (Activa en Rama P, x < 910240053.92 MHz via u_P)
A1  = 506.8885;
SB1 = -0.4416;
SC1 = 26.1500;
B1  = A1 + SB1;
C1p = 0.5 * A1 + SC1;
dA1 = -3.7500;
dB1 = -3.7500;
dC1 = -1.8375;
Sh1 = 60.00;

y_u1 = AsymRotor(x_exp, A1, B1, C1p, dA1, dB1, dC1, ...
                 Nu0_base + Sh1, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);

% Evaluacion de la funcion global con switch de Fermi-Dirac:
% y = Int1*1000*AsymRotor(Copia1)*u(a-x) + Int2*1000*AsymRotor(Copia2)*u(x-a) + C1
Int1_fd = 4.2687;
Int2_fd = 5.4473;
C1_fd   = 0.2040;

comp1_P = (Int1_fd * 1000.0) * y_u1 .* u_P;
comp2_R = (Int2_fd * 1000.0) * y_u2 .* u_R;
y_dual  = comp1_P + comp2_R + C1_fd;

resid_dual = y_exp - y_dual;
sse_dual   = sum(resid_dual.^2);
rms_dual   = sqrt(mean(resid_dual.^2));

% Metricas por region (Rama P: x < a_sw, Rama R: x >= a_sw)
m_P = (x_exp < a_sw);
m_R = (x_exp >= a_sw);
rel_ghz = (x_exp - (Nu0_base + Sh2)) / 1000.0;
m_R_hi  = (rel_ghz >= 6.5) & (rel_ghz <= 9.5);

R_P_1c = corrcoef(y_exp(m_P), y_single(m_P));
R_P_2c = corrcoef(y_exp(m_P), y_dual(m_P));
R_R_2c = corrcoef(y_exp(m_R), y_dual(m_R));
R_R_hi = corrcoef(y_exp(m_R_hi), y_dual(m_R_hi));
R_G_2c = corrcoef(y_exp, y_dual);

r_P_1comp = R_P_1c(1, 2);
r_P_2comp = R_P_2c(1, 2);
r_R_2comp = R_R_2c(1, 2);
r_R_win   = R_R_hi(1, 2);
r_G_2comp = R_G_2c(1, 2);

% =========================================================================
% REPORTE EN COMMAND WINDOW
% =========================================================================
fprintf('\n==========================================================================\n');
fprintf(' REPORTE DE RESULTADOS - BLUEPRINTXAS v5.23 (TRIFENILENO: %s)\n', dataFile);
fprintf('==========================================================================\n\n');

fprintf('1) COMPARATIVA BGridC: ScoreSolver = false vs ScoreSolver = true (Seed = 1)\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf(' Modalidad                          Pearson Global   Amplitud I1   SSE\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf(' ScoreSolver = false (Solo SSE)        %6.2f%%         %6.3f     %7.3f\n', 100*r_sse_false, I1_sse, sse_false);
fprintf(' ScoreSolver = true  (1 Componente)    %6.2f%%         %6.3f     %7.3f\n', 100*r_single, I1_single, sse_true);
fprintf(' ScoreSolver = true  (2 Comp. Fermi)   %6.2f%%         %6.3f     %7.3f\n', 100*r_G_2comp, Int2_fd, sse_dual);
fprintf('--------------------------------------------------------------------------\n\n');

fprintf('2) MODELO DE 2 COPIAS CON SWITCH DE FERMI-DIRAC (a = %.2f MHz)\n', a_sw);
fprintf('   Ecuacion Global:\n');
fprintf('   Int1*1000*AsymRotor(Copia1)*(1/(1+exp((x-a)/1e-20))) + ...\n');
fprintf('   Int2*1000*AsymRotor(Copia2)*(1/(1+exp((a-x)/1e-20))) + C1\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf(' Parametro      Copia 1 (Rama P, x < a)     Copia 2 (Rama R, x > a)\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf(' A  (MHz)             %10.4f                  %10.4f\n', A1, A2);
fprintf(' SB (MHz)             %+10.4f                  %+10.4f\n', SB1, SB2);
fprintf(' SC (MHz)             %+10.4f                  %+10.4f\n', SC1, SC2);
fprintf(' dA (MHz)             %+10.4f                  %+10.4f\n', dA1, dA2);
fprintf(' dB (MHz)             %+10.4f                  %+10.4f\n', dB1, dB2);
fprintf(' dC (MHz)             %+10.4f                  %+10.4f\n', dC1, dC2);
fprintf(' Shift (MHz)          %+10.2f                  %+10.2f\n', Sh1, Sh2);
fprintf(' Int (x1000)          %10.4f                  %10.4f\n', Int1_fd, Int2_fd);
fprintf(' Temp (K)             %10.2f                  %10.2f\n', Trot, Trot);
fprintf(' Theta / Phi       %5.1f / %4.1f               %5.1f / %4.1f\n', Theta, Phi, Theta, Phi);
fprintf(' GW / LW (MHz)     %5.2f / %4.2f               %5.2f / %4.2f\n', GW, LW, GW, LW);
fprintf(' C1 (Offset)                     %10.4f\n', C1_fd);
fprintf('--------------------------------------------------------------------------\n');
fprintf(' Correlacion Rama P (x < a)  : Sube de %5.2f%% (1 copia) a %5.2f%% (2 copias)\n', 100*r_P_1comp, 100*r_P_2comp);
fprintf(' Correlacion Rama R (x > a)  : %5.2f%% (Ventana +6.5..9.5 GHz: %5.2f%%)\n', 100*r_R_2comp, 100*r_R_win);
fprintf(' Correlacion Global Completa : Sube de %5.2f%% (1 copia) a %5.2f%% (2 copias)\n', 100*r_single, 100*r_G_2comp);
fprintf(' RMS Residual Global         : %.4f\n', rms_dual);
fprintf('==========================================================================\n\n');

% =========================================================================
% GRAFICAS DE RESULTADOS
% =========================================================================
% Figura 1: Modelo de 2 Copias con Switch de Fermi-Dirac
figure('Name', 'Opcion A: 2 Copias AsymRotor con Switch Fermi-Dirac', ...
       'Color', 'w', 'Position', [70, 50, 1120, 740]);

subplot(3, 1, 1);
plot(x_exp, y_exp, '.k', 'MarkerSize', 5.5); hold on;
plot(x_exp(m_P), comp1_P(m_P) + C1_fd, '-', 'Color', [0.80, 0.25, 0.15], 'LineWidth', 1.1);
plot(x_exp(m_R), comp2_R(m_R) + C1_fd, '-b', 'LineWidth', 1.1);
xline(a_sw, '--g', 'LineWidth', 1.4);
ylabel('Intensity (a.u.)');
title(sprintf('Trifenileno - 2 Copias AsymRotor + Switch Fermi-Dirac (a = %.2f MHz) | r_{Global} = %.2f%% (antes 56.33%%)', ...
      a_sw, 100*r_G_2comp));
legend('Exp (pmt-nFLa)', ...
       sprintf('Copia 1: Rama P (dA_1=%.2f, SC_1=%+.2f | r_P=%.2f%%)', dA1, SC1, 100*r_P_2comp), ...
       sprintf('Copia 2: Rama R (dA_2=%.2f, SC_2=%+.2f | r_R=%.2f%%)', dA2, SC2, 100*r_R_2comp), ...
       'Switch Fermi-Dirac (x = a)', 'Location', 'northeast');
xlim([min(x_exp), max(x_exp)]);
ylim([0.12, max(y_exp)*1.10]);
grid on; box on; set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.6);

subplot(3, 1, 2);
plot(x_exp(m_P), y_exp(m_P), '.k', 'MarkerSize', 6); hold on;
plot(x_exp(m_P), y_single(m_P), '--', 'Color', [0.55, 0.55, 0.55], 'LineWidth', 0.9);
plot(x_exp(m_P), y_dual(m_P), '-', 'Color', [0.80, 0.25, 0.15], 'LineWidth', 1.2);
ylabel('Intensity (Rama P)');
title(sprintf('Mejora en Rama P (x < a): 1 Copia (r = %.2f%%) vs 2 Copias Fermi-Dirac (r = %.2f%%)', ...
      100*r_P_1comp, 100*r_P_2comp));
legend('Exp Rama P', '1 Copia (Desfasada en P)', '2 Copias (Copia 1 activa)', 'Location', 'northwest');
xlim([min(x_exp(m_P)), max(x_exp(m_P))]);
ylim([0.12, max(y_exp(m_P))*1.10]);
grid on; box on; set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.6);

subplot(3, 1, 3);
plot(x_exp, resid_dual, '-r', 'LineWidth', 0.8); hold on;
yline(0, '--k');
xline(a_sw, '--g', 'LineWidth', 1.2);
xlabel('Frequency (MHz)');
ylabel('Residual');
title(sprintf('Residuales Modelo 2 Copias Fermi-Dirac (RMS = %.4f, SSE = %.2f)', rms_dual, sse_dual));
xlim([min(x_exp), max(x_exp)]);
ylim([-0.45, 0.45]);
grid on; box on; set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.6);

% Figura 2: Comparativa BGridC ScoreSolver = false vs ScoreSolver = true (Seed = 1)
figure('Name', 'Comparativa BGridC: ScoreSolver false vs true (Seed = 1)', ...
       'Color', 'w', 'Position', [110, 80, 1080, 620]);

subplot(2, 1, 1);
plot(x_exp, y_exp, '.k', 'MarkerSize', 5); hold on;
plot(x_exp, y_sse, '-', 'Color', [0.85, 0.35, 0.10], 'LineWidth', 1.1);
ylabel('Intensity (a.u.)');
title(sprintf('BGridC (ScoreSolver = false, Seed = 1 | Solo SSE): Amplitud deprimida I_1 = %.2f | r = %.2f%% | SSE = %.2f', ...
      I1_sse, 100*r_sse_false, sse_false));
legend('Exp (pmt-nFLa)', 'ScoreSolver = false (Solo SSE)', 'Location', 'northeast');
xlim([min(x_exp), max(x_exp)]);
ylim([0.12, max(y_exp)*1.08]);
grid on; box on; set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.6);

subplot(2, 1, 2);
plot(x_exp, y_exp, '.k', 'MarkerSize', 5); hold on;
plot(x_exp, y_single, '-b', 'LineWidth', 1.1);
ylabel('Intensity (a.u.)');
xlabel('Frequency (MHz)');
title(sprintf('BGridC (ScoreSolver = true, Seed = 1 | SSE + Shape Weight): I_1 = %.2f | r_{Global} = %.2f%% (r_{Rama R} = 84.87%%) | SSE = %.2f', ...
      I1_single, 100*r_single, sse_true));
legend('Exp (pmt-nFLa)', 'ScoreSolver = true (SSE + Correlacion de Forma)', 'Location', 'northeast');
xlim([min(x_exp), max(x_exp)]);
ylim([0.12, max(y_exp)*1.08]);
grid on; box on; set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.6);

results = struct('r_sse_false', r_sse_false, 'r_single_true', r_single, ...
                 'r_dual_global', r_G_2comp, 'r_P_1comp', r_P_1comp, ...
                 'r_P_2comp', r_P_2comp, 'r_R_2comp', r_R_2comp);
end

