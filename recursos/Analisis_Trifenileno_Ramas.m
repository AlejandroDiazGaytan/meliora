%% Analisis_Trifenileno_Ramas.m
% =========================================================================
% ANALISIS ESPECTROSCOPICO Y DESCOMPOSICION POR RAMAS: TRIFENILENO (C18H12)
% Comparacion del espectro experimental LIF (pmt-nFLa) vs simulacion AsymRotor
%
% Division de Ciencias e Ingenierias (Campus Leon) - Universidad de Guanajuato
% Autor: Ing. Miguel Alejandro Diaz Gaytan (Alex DG)
% Directores de Tesis: Dr. Leonardo Alvarez-Valtierra y Dr. Mario U. Delgado-Jaime
% Colaboradores: Dr. Jaime Valdez-Ruvalcaba, Dr. Alejandro Gil-Villegas Montiel
% =========================================================================
%
% Mis notas para la discusion de resultados con el posdoc:
% -------------------------------------------------------------------------
% 1. El trifenileno es un rotor simetrico oblato plano de simetria D3h (kappa -> +1).
% 2. En el modelo plano rigido estricto:
%      C = (A * B) / (A + B) = B / 2
%    Esto anula el factor 2*(B - 2*C) = 0, provocando el colapso de los
%    desdoblamientos rotacionales de los bloques de Wang. Por eso un ajuste
%    plano estricto solo logra r ~ 13.45% a 16.77%.
% 3. La banda observada a ~329.35 nm (nu0 = 910,239.9 GHz / 30,362.34 cm-1)
%    es un falso origen vibronico Herzberg-Teller inducido por un modo degenerado e'.
% 4. El momento angular vibronic (l = +/- 1) interactua con la rotacion molecular
%    generando acoplamiento de Coriolis de primer orden, el cual se absorbe en una
%    constante rotacional axial efectiva:
%      C_eff = C_plano * (1 - zeta_eff)  -->  C_eff = A/2 + SC
%    con SC = +26.40 MHz, lo que da C_eff = 279.84 MHz frente a C_plano = 253.33 MHz.
% 5. Al incorporar SC = +26.40 MHz:
%    - Se reproducen fielmente los sub-picos periodicos de 90.4 MHz (~0.10 cm-1).
%    - La correlacion en la Rama R salta a r = 84.64% en toda la rama (+2 a +10 GHz)
%      y alcanza r = 92.90% (Coseno = 95.57%) en la ventana de alta resolucion (+6.5 a +9.5 GHz).
% 6. La Rama P presenta una correlacion menor (r ~ 28.8%) debido a la modulacion
%    de intensidad inducida por interferencia vibronica tipo Herman-Wallis.
% 7. La region central presenta el hueco caracteristico de transicion Tipo B (E' <- A1').
% =========================================================================

clear; clc; close all;

fprintf('========================================================================\n');
fprintf('  ANALISIS ESPECTROSCOPICO DE TRIFENILENO (C18H12) - MODO ASYMROTOR\n');
fprintf('  Datos experimentales: pmt-nFLa (LIF chorro supersonico, DCI-UG)\n');
fprintf('  Modelo: Rotor cuasi-oblato D3h con perturbacion de Coriolis efectiva\n');
fprintf('========================================================================\n\n');

%% 1. Verificacion de rutas y carga de datos experimentales
% Asegurar que la funcion AsymRotor y los datos esten accesibles en el path
if ~exist('AsymRotor', 'file')
    if exist(fullfile(pwd, 'recursos'), 'dir')
        addpath(fullfile(pwd, 'recursos'));
    elseif exist(fullfile(fileparts(mfilename('fullpath')), 'recursos'), 'dir')
        addpath(fullfile(fileparts(mfilename('fullpath')), 'recursos'));
    end
end

dataFile = 'pmt-nFLa';
if ~exist(dataFile, 'file')
    if exist(fullfile('recursos', 'pmt-nFLa'), 'file')
        dataFile = fullfile('recursos', 'pmt-nFLa');
    elseif exist(fullfile('..', 'recursos', 'pmt-nFLa'), 'file')
        dataFile = fullfile('..', 'recursos', 'pmt-nFLa');
    else
        error('No se encontro el archivo experimental %s en la carpeta de trabajo.', dataFile);
    end
end

fprintf('Cargando archivo experimental: %s ...\n', dataFile);
raw = load(dataFile);
x_exp = raw(:, 1);  % Frecuencia experimental en MHz
y_exp = raw(:, 2);  % Intensidad experimental cruda (LIF)
nPts  = length(x_exp);

% Constantes de conversion
C_LIGHT_CM = 29979.2458; % MHz por cm-1
x_cm = x_exp / C_LIGHT_CM;

fprintf('  Puntos espectrales registrados: %d\n', nPts);
fprintf('  Rango espectral: [%.2f, %.2f] MHz  |  [%.3f, %.3f] cm-1\n\n', ...
        min(x_exp), max(x_exp), min(x_cm), max(x_cm));

%% 2. Parametros rotacionales y fisicos calibrados
% Estado Fundamental S0
A0      = 506.88854;                % MHz
SB      = -0.44159;                 % Desviacion de simetria B = A + SB
SC      = +26.39768;                % Acoplamiento de Coriolis efectivo C = A/2 + SC
B0      = A0 + SB;                  % 506.44695 MHz
C0_eff  = (A0 / 2.0) + SC;          % 279.84195 MHz (Coriolis efectiva)
C0_plan = (A0 * B0) / (A0 + B0);    % 253.33230 MHz (Plano rigido de referencia)

% Cambios de Excitacion S1 <- S0
dA      = -3.09507;                 % MHz
dB      = -3.09941;                 % MHz
dC      = -1.46120;                 % MHz

% Origen vibronico y condiciones del haz supersonico
Nu0_ref = 910238257.8;              % MHz (origen de referencia en BlueprintXAS)
Shift   = +1657.39;                 % MHz (desplazamiento de calibracion)
Nu0     = Nu0_ref + Shift;          % 910,239,915.19 MHz (30,362.335 cm-1)
Nu0_cm  = Nu0 / C_LIGHT_CM;

Trot     = 3.04;                    % K (temperatura rotacional principal)
Theta    = 44.8;                    % Grados (banda hibrida en plano aromatico e': mu_a^2 = mu_b^2 = 0.5)
Phi      = 0.0;                     % Grados
GW       = 31.53;                   % MHz (FWHM Gauss Doppler / instrumental)
LW       = 3.10;                    % MHz (FWHM Lorentz tiempo de vida S1)
JMax     = 50;                      % Nivel rotacional maximo
Unit     = 1;                       % 1 = MHz
SpinAxis = 2;                       % Simetria eje b (C2)
SpinOdd  = 2.0;                     % Peso estadistico nuclear
T2       = 12.0;                    % K (cola termica secundaria del jet supersonico)
W2       = 0.0501;                  % Fraccion de la componente T2

kappa = (2*B0 - A0 - C0_eff) / (A0 - C0_eff);

%% 3. Evaluacion del modelo de AsymRotor y escalamiento lineal
fprintf('Calculando simulacion cuantica mediante AsymRotor...\n');
[y_unit, Lines] = AsymRotor(x_exp, A0, B0, C0_eff, dA, dB, dC, ...
                            Nu0, Trot, Theta, Phi, GW, LW, ...
                            JMax, Unit, SpinAxis, SpinOdd, T2, W2);

% Ajuste lineal de linea base e intensidad: y_fit = I1 * y_unit + Offset
M_lin = [y_unit, ones(nPts, 1)];
coef  = M_lin \ y_exp;
Int_scale = coef(1);
Offset    = coef(2);
I1_ref    = Int_scale / 1000.0;     % En convencion de ecuacion BlueprintXAS

y_fit     = Int_scale * y_unit + Offset;
residuales = y_exp - y_fit;
rms_fit   = sqrt(mean(residuales.^2));

%% 4. Evaluacion de correlaciones por regiones y ramas espectrales
% Frecuencia relativa al origen vibronico en GHz
rel_ghz = (x_exp - Nu0) / 1000.0;

% Mascaras de las ventanas definidas en el poster
mask_Global = true(nPts, 1);
mask_P_post = (rel_ghz >= -9.5) & (rel_ghz <= -6.5);
mask_P_wide = (rel_ghz >= -10.0) & (rel_ghz <= -2.0);
mask_Q_cent = (rel_ghz >= -1.5) & (rel_ghz <= 1.5);
mask_R_post = (rel_ghz >= 6.5) & (rel_ghz <= 9.5);
mask_R_wide = (rel_ghz >= 2.0) & (rel_ghz <= 10.0);

% Calculo de coeficientes de Pearson (r)
r_global = calc_r(y_exp(mask_Global), y_fit(mask_Global));
r_P_post = calc_r(y_exp(mask_P_post), y_fit(mask_P_post));
r_P_wide = calc_r(y_exp(mask_P_wide), y_fit(mask_P_wide));
r_Q_cent = calc_r(y_exp(mask_Q_cent), y_fit(mask_Q_cent));
r_R_post = calc_r(y_exp(mask_R_post), y_fit(mask_R_post));
r_R_wide = calc_r(y_exp(mask_R_wide), y_fit(mask_R_wide));

% Calculo de similitud coseno (sustrayendo linea base estimada)
cos_global = calc_cos(y_exp(mask_Global), y_fit(mask_Global), Offset);
cos_P_post = calc_cos(y_exp(mask_P_post), y_fit(mask_P_post), Offset);
cos_P_wide = calc_cos(y_exp(mask_P_wide), y_fit(mask_P_wide), Offset);
cos_Q_cent = calc_cos(y_exp(mask_Q_cent), y_fit(mask_Q_cent), Offset);
cos_R_post = calc_cos(y_exp(mask_R_post), y_fit(mask_R_post), Offset);
cos_R_wide = calc_cos(y_exp(mask_R_wide), y_fit(mask_R_wide), Offset);

%% 5. Impresion limpia de resultados en Command Window
fprintf('\n========================================================================\n');
fprintf('  PARAMETROS DEL AJUSTE (FORMATO BLUEPRINTXAS / ASYMROTOR)\n');
fprintf('========================================================================\n');
fprintf('  Ecuacion de Ajuste:\n');
fprintf('    y(x) = I1*1000*AsymRotor(x, A, A+SB, A/2+SC, dA, dB, dC, Nu0+Shift, ...) + C1\n\n');
fprintf('  Constantes Estado Fundamental S0:\n');
fprintf('    A"           = %10.4f MHz\n', A0);
fprintf('    B"           = %10.4f MHz   (SB = %+8.4f MHz)\n', B0, SB);
fprintf('    C" (C_eff)   = %10.4f MHz   (SC = %+8.4f MHz, Coriolis D3h)\n', C0_eff, SC);
fprintf('    C_plano ref  = %10.4f MHz   (A*B/(A+B))\n', C0_plan);
fprintf('    Ray kappa    = %+10.4f     (Rotor simetrico oblato casi puro)\n\n', kappa);
fprintf('  Cambios de Excitacion S1 <- S0:\n');
fprintf('    dA           = %+10.4f MHz   -->  A'' = %.4f MHz\n', dA, A0 + dA);
fprintf('    dB           = %+10.4f MHz   -->  B'' = %.4f MHz\n', dB, B0 + dB);
fprintf('    dC           = %+10.4f MHz   -->  C'' = %.4f MHz\n\n', dC, C0_eff + dC);
fprintf('  Condiciones del Haz y Espectrometro:\n');
fprintf('    Nu0 absoluto = %12.2f MHz (%.3f cm-1 | Shift = %+.2f MHz)\n', Nu0, Nu0_cm, Shift);
fprintf('    Trot         = %10.2f K     (T2 = %.1f K, W2 = %.4f)\n', Trot, T2, W2);
fprintf('    GW / LW      =   %6.2f / %.2f MHz\n', GW, LW);
fprintf('    Dipolo       =   Theta = %.1f deg, Phi = %.1f deg (Banda hibrida en plano e'')\n', Theta, Phi);
fprintf('    I1 / C1      =   %6.3f / %.4f  (RMS residual = %.4f)\n', I1_ref, Offset, rms_fit);
fprintf('------------------------------------------------------------------------\n');
fprintf('  METRICAS DE CORRELACION ESPECTROSCOPICA POR RAMAS:\n');
fprintf('  Region / Rama                Ventana Relativa Nu0    Pearson (r)    Similitud Coseno\n');
fprintf('  ----------------------------------------------------------------------\n');
fprintf('  Espectro Global Completo     [-18.5 a +18.5 GHz]       %5.2f %%           %5.2f %%\n', r_global, cos_global);
fprintf('  Rama R (Ventana Poster)      [ +6.5 a  +9.5 GHz]       %5.2f %%           %5.2f %%\n', r_R_post, cos_R_post);
fprintf('  Rama R Amplia                [ +2.0 a +10.0 GHz]       %5.2f %%           %5.2f %%\n', r_R_wide, cos_R_wide);
fprintf('  Region Central (Rama Q)      [ -1.5 a  +1.5 GHz]       %5.2f %%           %5.2f %%\n', r_Q_cent, cos_Q_cent);
fprintf('  Rama P (Ventana Poster)      [ -9.5 a  -6.5 GHz]       %5.2f %%           %5.2f %%\n', r_P_post, cos_P_post);
fprintf('  Rama P Amplia                [-10.0 a  -2.0 GHz]       %5.2f %%           %5.2f %%\n', r_P_wide, cos_P_wide);
fprintf('========================================================================\n\n');

%% 6. Figura 1: Espectro Completo Experimental vs Simulacion (Layout General)
fig1 = figure('Name', 'Trifenileno - Comparacion Espectro Completo (pmt-nFLa vs AsymRotor)', ...
              'Color', 'w', 'Position', [80, 80, 1100, 720]);

% Panel Superior: Comparacion de Espectros
subplot(3, 1, [1, 2]);
plot(x_cm, y_exp, '-k', 'LineWidth', 0.9, 'DisplayName', 'Experimental (pmt-nFLa)');
hold on;
plot(x_cm, y_fit, '-b', 'LineWidth', 1.1, ...
     'DisplayName', sprintf('Simulacion AsymRotor (r = %.2f%%, Coseno = %.2f%%)', r_global, cos_global));

% Anotaciones de ramas
xline(Nu0_cm, '--g', 'LineWidth', 1.2, 'DisplayName', sprintf('Origen Vibronico \\nu_0 = %.2f cm^{-1}', Nu0_cm));
text(Nu0_cm, max(y_exp)*0.95, '  \nu_0 (30,362.34 cm^{-1})', 'Color', [0, 0.5, 0], 'FontSize', 9, 'FontWeight', 'bold');
text(30362.05, max(y_exp)*0.85, 'Rama P (\DeltaJ = -1)', 'Color', [0.5, 0, 0.5], 'FontSize', 9, 'FontWeight', 'bold');
text(30362.62, max(y_exp)*0.85, 'Rama R (\DeltaJ = +1)', 'Color', [0, 0.35, 0.8], 'FontSize', 9, 'FontWeight', 'bold');

ylabel('Intensidad de Fluorescencia LIF (u.a.)', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Trifenileno (C_{18}H_{12}) · Transicion Vibronica a 329.35 nm (S_1 \\leftarrow S_0)\nA" = %.3f MHz, B" = %.3f MHz, C"_{eff} = %.3f MHz (SC = %+5.2f MHz) | T_{rot} = %.2f K', ...
      A0, B0, C0_eff, SC, Trot), 'FontSize', 11, 'FontWeight', 'bold');
xlim([min(x_cm), max(x_cm)]);
ylim([0.15, max(y_exp) * 1.10]);
legend('Location', 'northeast', 'FontSize', 9);
grid on; box on;

% Panel Inferior: Residuales (Exp - Sim)
subplot(3, 1, 3);
plot(x_cm, residuales, '-r', 'LineWidth', 0.8, 'DisplayName', sprintf('Residuales (RMS = %.4f)', rms_fit));
hold on;
yline(0, '--k', 'LineWidth', 0.8);
xlabel('Frecuencia (cm^{-1})', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Residual (u.a.)', 'FontSize', 10, 'FontWeight', 'bold');
xlim([min(x_cm), max(x_cm)]);
ylim([-0.40, 0.40]);
legend('Location', 'northeast', 'FontSize', 9);
grid on; box on;

%% 7. Figura 2: Descomposicion por Ramas y Detalle de Sub-picos (Layout Poster UG)
fig2 = figure('Name', 'Trifenileno - Analisis por Ramas y Resolucion de Sub-picos (Poster UG)', ...
              'Color', 'w', 'Position', [120, 60, 1200, 780]);

% Ventanas en cm-1
x_P_cm = x_cm(mask_P_post);
y_exp_P = y_exp(mask_P_post);
y_fit_P = y_fit(mask_P_post);

x_Q_cm = x_cm(mask_Q_cent);
y_exp_Q = y_exp(mask_Q_cent);
y_fit_Q = y_fit(mask_Q_cent);

x_R_cm = x_cm(mask_R_post);
y_exp_R = y_exp(mask_R_post);
y_fit_R = y_fit(mask_R_post);

% Subplot (2, 3, 1): Rama P - Espectro Experimental Crudo
subplot(2, 3, 1);
plot(x_P_cm, y_exp_P, '-k', 'LineWidth', 1.0);
title(sprintf('Rama P (-9.5 a -6.5 GHz)\nInterferencia Vibronica Herman-Wallis'), 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Experimental Crudo', 'FontWeight', 'bold');
xlabel('Frecuencia (cm^{-1})');
xlim([min(x_P_cm), max(x_P_cm)]);
ylim([0.18, max(y_exp_P)*1.10]);
grid on; box on;

% Subplot (2, 3, 2): Region Central Q - Espectro Experimental Crudo
subplot(2, 3, 2);
plot(x_Q_cm, y_exp_Q, '-k', 'LineWidth', 1.0);
title(sprintf('Region Central \\nu_0 (-1.5 a +1.5 GHz)\nHueco Central Tipo B (E'' \\leftarrow A_1'')'), 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Experimental Crudo', 'FontWeight', 'bold');
xlabel('Frecuencia (cm^{-1})');
xlim([min(x_Q_cm), max(x_Q_cm)]);
ylim([0.18, max(y_exp_Q)*1.10]);
grid on; box on;

% Subplot (2, 3, 3): Rama R - Espectro Experimental Crudo
subplot(2, 3, 3);
plot(x_R_cm, y_exp_R, '-k', 'LineWidth', 1.0);
title(sprintf('Rama R (+6.5 a +9.5 GHz)\nr = %.2f%% (Coseno = %.2f%%)', r_R_post, cos_R_post), ...
      'FontSize', 10, 'FontWeight', 'bold', 'Color', [0, 0.35, 0.8]);
ylabel('Experimental Crudo', 'FontWeight', 'bold');
xlabel('Frecuencia (cm^{-1})');
xlim([min(x_R_cm), max(x_R_cm)]);
ylim([0.18, max(y_exp_R)*1.10]);
grid on; box on;

% Subplot (2, 3, 4): Rama P - Comparacion y Sub-picos
subplot(2, 3, 4);
plot(x_P_cm, y_exp_P, '-k', 'LineWidth', 0.9, 'DisplayName', 'Exp');
hold on;
plot(x_P_cm, y_fit_P, '-b', 'LineWidth', 1.1, 'DisplayName', sprintf('AsymRotor (r = %.2f%%)', r_P_post));
title('Resolucion de Sub-picos (~90.4 MHz)', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('Exp vs AsymRotor', 'FontWeight', 'bold');
xlabel('Frecuencia (cm^{-1})');
xlim([min(x_P_cm), max(x_P_cm)]);
ylim([0.18, max(y_exp_P)*1.10]);
legend('Location', 'northeast', 'FontSize', 8);
grid on; box on;

% Subplot (2, 3, 5): Region Central Q - Comparacion y Ancho Instrumental
subplot(2, 3, 5);
plot(x_Q_cm, y_exp_Q, '-k', 'LineWidth', 0.9, 'DisplayName', 'Exp');
hold on;
plot(x_Q_cm, y_fit_Q, '-b', 'LineWidth', 1.1, 'DisplayName', sprintf('AsymRotor (r = %.2f%%)', r_Q_cent));
title('Perfil Central (~37.2 MHz GW)', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('Exp vs AsymRotor', 'FontWeight', 'bold');
xlabel('Frecuencia (cm^{-1})');
xlim([min(x_Q_cm), max(x_Q_cm)]);
ylim([0.18, max(y_exp_Q)*1.10]);
legend('Location', 'northeast', 'FontSize', 8);
grid on; box on;

% Subplot (2, 3, 6): Rama R - Comparacion y Sub-picos
subplot(2, 3, 6);
plot(x_R_cm, y_exp_R, '-k', 'LineWidth', 0.9, 'DisplayName', 'Exp');
hold on;
plot(x_R_cm, y_fit_R, '-b', 'LineWidth', 1.2, 'DisplayName', sprintf('AsymRotor (r = %.2f%%)', r_R_post));
title('Resolucion de Sub-picos (~90.4 MHz)', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('Exp vs AsymRotor', 'FontWeight', 'bold');
xlabel('Frecuencia (cm^{-1})');
xlim([min(x_R_cm), max(x_R_cm)]);
ylim([0.18, max(y_exp_R)*1.10]);
legend('Location', 'northeast', 'FontSize', 8);
grid on; box on;

fprintf('Graficas generadas exitosamente (Figura 1 y Figura 2).\n');

%% Funciones auxiliares locales
function r = calc_r(y1, y2)
    C = corrcoef(y1, y2);
    r = C(1, 2) * 100.0;
    if isnan(r), r = 0.0; end
end

function c = calc_cos(y1, y2, offset)
    ye0 = max(0, y1 - offset);
    ys0 = max(0, y2 - offset);
    norm_prod = norm(ye0) * norm(ys0);
    if norm_prod > 1e-12
        c = (dot(ye0, ys0) / norm_prod) * 100.0;
    else
        c = 0.0;
    end
end
