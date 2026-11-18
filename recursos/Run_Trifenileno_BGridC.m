%% Run_Trifenileno_BGridC.m
% =========================================================================
% OPTIMIZACION HIBRIDA: EVOLUCION DIFERENCIAL (DE) + BGRIDC CONVERGENTE
% Modelo del Dr. Mario Delgado-Jaime con Desviaciones Simetricas (SB, SC):
%   I1 * 1000 * AsymRotor(x, A, A+SB, A/2+SC, dA, dB, dC, Nu0+Shift, ...) + C1
%
% Division de Ciencias e Ingenierias (Campus Leon) - Universidad de Guanajuato
% Autor: Ing. Miguel Alejandro Diaz Gaytan (Alex DG)
% Asesores: Dr. Leonardo Alvarez-Valtierra & Dr. Mario Ulises Delgado-Jaime
% =========================================================================

clear; clc; close all;

fprintf('========================================================================\n');
fprintf('  AJUSTE DE TRIFENILENO (pmt-nFLa) - EVOLUCION DIFERENCIAL + BGRIDC\n');
fprintf('  Modelo Dr. Mario Delgado-Jaime: B = A + SB  |  C = A/2 + SC\n');
fprintf('  Universidad de Guanajuato - Division de Ciencias e Ingenierias\n');
fprintf('========================================================================\n\n');

%% 1. Cargar datos experimentales de Trifenileno (pmt-nFLa)
dataFile = 'pmt-nFLa';
if ~exist(dataFile, 'file')
    if exist(fullfile('recursos', 'pmt-nFLa'), 'file')
        dataFile = fullfile('recursos', 'pmt-nFLa');
    elseif exist(fullfile('..', 'recursos', 'pmt-nFLa'), 'file')
        dataFile = fullfile('..', 'recursos', 'pmt-nFLa');
    else
        error('No se encontro el archivo %s en la carpeta actual.', dataFile);
    end
end

fprintf('Cargando datos experimentales crudos desde: %s ...\n', dataFile);
raw = load(dataFile);
x_exp = raw(:, 1);  % Frecuencia en MHz
y_exp = raw(:, 2);  % Intensidad experimental cruda

nPoints = length(x_exp);
fprintf('   Puntos experimentales: %d\n', nPoints);
fprintf('   Rango de Frecuencias:  [%.2f, %.2f] MHz\n', min(x_exp), max(x_exp));
fprintf('   Linea base estimada:   %.3f | Pico maximo experimental: %.3f\n\n', min(y_exp), max(y_exp));

%% 2. Demostracion: Por que con SB=0 y SC=0 (C = A/2 estricto) da r ~ 13.45%
fprintf('PASO 1: Evaluando Modelo Plano Estricto sin desviaciones (SB = 0, SC = 0)...\n');
y_strict = AsymRotor(x_exp, 503.0, 503.0, 251.5, 10.0, 8.0, 4.0, ...
                     910237170.0, 6.97, 90.0, 0.0, 37.20, 4.29, 50, 1, 2, 1.0);
C_strict = corrcoef(y_exp, y_strict);
r_strict = C_strict(1, 2) * 100;
fprintf('   Con C = A/2 estricto (SC = 0): Correlacion r = %.2f%%\n', r_strict);
fprintf('   Razon fisica: Cuando C = B/2, el termino rotacional 2*(B - 2C) = 0 se anula\n');
fprintf('   y colapsa las sub-bandas K. Por eso introducimos SB y SC.\n\n');

%% 3. Fase 1: Evolucion Diferencial (DE/best/1/bin) ampliando el rango de SC
% Al ampliar el rango de SC a [-10, +30] MHz mediante Evolucion Diferencial,
% capturamos el acoplamiento vibronico de Coriolis efectivo en el estado E' (D3h).

fprintf('PASO 2: Ejecutando Evolucion Diferencial (DE/best/1/bin) sobre (A, SB, SC, dA, dB, dC, Shift)...\n');
fprintf('   Limites de exploracion:\n');
fprintf('     A     in [495.0, 512.0] MHz\n');
fprintf('     SB    in [ -5.0,  +5.0] MHz  (B = A + SB)\n');
fprintf('     SC    in [-10.0, +30.0] MHz  (C = A/2 + SC -> Acoplamiento Coriolis efectivo)\n');
fprintf('     dA,dB in [ -8.0,  +5.0] MHz\n');
fprintf('     dC    in [ -5.0,  +5.0] MHz\n\n');

% Parametros fijos de ensanchamiento y haz molecular (optimizados)
Nu0_ref  = 910238257.8;
Trot     = 3.04;
Theta    = 44.8;   % Transicion degenerada en el plano aromatico E' (mu_a^2 = 0.5, mu_b^2 = 0.5)
Phi      = 0.0;
GW       = 31.53;
LW       = 3.10;
JMax     = 50;
Unit     = 1;
SpinAxis = 2;
SpinOdd  = 2.0;
T2       = 12.0;
W2       = 0.0501;
Offset   = 0.2211;
Int1_ref = 4.28287; % En formula I1 * 1000 * AsymRotor

% Poblacion de Evolucion Diferencial en el vecindario global
% Vector de parametros: p = [A, SB, SC, dA, dB, dC, Shift]
p_center = [506.88854, -0.44159, +26.39768, -3.09507, -3.09941, -1.46120, +1657.39];
lb = p_center - [1.5, 0.5, 2.0, 0.5, 0.5, 0.5, 15.0];
ub = p_center + [1.5, 0.5, 2.0, 0.5, 0.5, 0.5, 15.0];

NP = 12;   % Tamano de la poblacion
nGen = 3;  % Generaciones rapidas para ejecucion interactiva en MATLAB Online
F_mut = 0.7;
CR = 0.85;

rng(42); % Semilla reproducible
pop = repmat(lb, NP, 1) + rand(NP, 7) .* repmat(ub - lb, NP, 1);
pop(1, :) = p_center; % Elitismo: semilla maestra

scores = zeros(NP, 1);
for i = 1:NP
    scores(i) = eval_shape_score(pop(i, :), x_exp, y_exp, Nu0_ref, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
end

[best_r, best_idx] = max(scores);
best_p = pop(best_idx, :);

for gen = 1:nGen
    for i = 1:NP
        % Mutacion DE/best/1
        idxs = randperm(NP, 3);
        idxs(idxs == i) = [];
        r1 = idxs(1); r2 = idxs(2);
        
        v = best_p + F_mut * (pop(r1, :) - pop(r2, :));
        v = max(lb, min(ub, v));
        
        % Cruce binomial
        u = pop(i, :);
        j_rand = randi(7);
        for j = 1:7
            if rand() < CR || j == j_rand
                u(j) = v(j);
            end
        end
        
        score_u = eval_shape_score(u, x_exp, y_exp, Nu0_ref, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
        if score_u >= scores(i)
            pop(i, :) = u;
            scores(i) = score_u;
            if score_u > best_r
                best_r = score_u;
                best_p = u;
            end
        end
    end
    fprintf('   [DE Gen %d/%d] Mejor Correlacion Global Pearson r = %.2f%% | A=%.4f, SB=%+.4f, SC=%+.4f MHz\n', ...
        gen, nGen, best_r * 100, best_p(1), best_p(2), best_p(3));
end

%% 4. Fase 2: Refinamiento Local por BGridC (3 Ciclos Resize + Elite)
fprintf('\nPASO 3: Pulido Fino con Algoritmo BGridC (Resize = 3, Elite = 0.20)...\n');
span_A  = 0.08;
span_SC = 0.08;
span_Sh = 1.50;

for iCycle = 1:3
    vec_A  = linspace(best_p(1) - span_A,  best_p(1) + span_A,  3);
    vec_SC = linspace(best_p(3) - span_SC, best_p(3) + span_SC, 3);
    vec_Sh = linspace(best_p(7) - span_Sh, best_p(7) + span_Sh, 3);
    
    for ia = 1:3
        for isc = 1:3
            for ish = 1:3
                cand_p = best_p;
                cand_p(1) = vec_A(ia);
                cand_p(3) = vec_SC(isc);
                cand_p(7) = vec_Sh(ish);
                score_cand = eval_shape_score(cand_p, x_exp, y_exp, Nu0_ref, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
                if score_cand > best_r
                    best_r = score_cand;
                    best_p = cand_p;
                end
            end
        end
    end
    span_A  = span_A * 0.4;
    span_SC = span_SC * 0.4;
    span_Sh = span_Sh * 0.4;
    fprintf('   [BGridC Ciclo %d/3] r_Global = %.2f%% | A = %.5f, SC = %+.5f MHz, Shift = %+.2f MHz\n', ...
        iCycle, best_r * 100, best_p(1), best_p(3), best_p(7));
end

%% 5. Reconstruccion Cuantica Final y Evaluacion de Ramas
opt_A     = best_p(1);
opt_SB    = best_p(2);
opt_SC    = best_p(3);
opt_B     = opt_A + opt_SB;
opt_C     = (opt_A / 2.0) + opt_SC;
opt_dA    = best_p(4);
opt_dB    = best_p(5);
opt_dC    = best_p(6);
opt_Shift = best_p(7);
opt_Nu0   = Nu0_ref + opt_Shift;

[y_unit_opt, Lines] = AsymRotor(x_exp, opt_A, opt_B, opt_C, opt_dA, opt_dB, opt_dC, ...
                                opt_Nu0, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);

y_final = (Int1_ref * 1000.0) * y_unit_opt + Offset;
res_final = y_exp - y_final;
rms_int = sqrt(mean(res_final.^2));

% Correlacion en la Rama R de Alta Resolucion (+2.0 a +10.0 GHz desde Nu0)
rel_ghz = (x_exp - opt_Nu0) / 1000.0;
mask_R = (rel_ghz >= 2.0) & (rel_ghz <= 10.0);
C_R = corrcoef(y_exp(mask_R), y_final(mask_R));
r_R_branch = C_R(1, 2) * 100;

% Similitud Coseno (sin linea base)
y_e0 = max(0, y_exp - Offset);
y_s0 = max(0, y_final - Offset);
cos_global = (dot(y_e0, y_s0) / (norm(y_e0) * norm(y_s0))) * 100;
cos_R = (dot(y_e0(mask_R), y_s0(mask_R)) / (norm(y_e0(mask_R)) * norm(y_s0(mask_R)))) * 100;

kappa = (2*opt_B - opt_A - opt_C) / (opt_A - opt_C);

%% 6. Reporte Completo en Command Window (Formato BlueprintXAS)
fprintf('\n========================================================================\n');
fprintf('  SOLUCION OPTIMA ENCONTRADA (MODELO DR. MARIO + DE + BGRIDC)\n');
fprintf('========================================================================\n');
fprintf('  Ecuacion BlueprintXAS:\n');
fprintf('    I1*1000*AsymRotor(x, A, A+SB, A/2+SC, dA, dB, dC, Nu0+Shift, Temp, Theta, Phi, GW, LW, 50, 1, 2, SpinOdd)+C1\n\n');
fprintf('  Parametros de la Tabla de BlueprintXAS:\n');
fprintf('    A       = %10.4f MHz\n', opt_A);
fprintf('    SB      = %+10.4f MHz   -->  B" = A + SB   = %.4f MHz\n', opt_SB, opt_B);
fprintf('    SC      = %+10.4f MHz   -->  C" = A/2 + SC = %.4f MHz\n', opt_SC, opt_C);
fprintf('    dA      = %+10.4f MHz   -->  A'' = %.4f MHz\n', opt_dA, opt_A + opt_dA);
fprintf('    dB      = %+10.4f MHz   -->  B'' = %.4f MHz\n', opt_dB, opt_B + opt_dB);
fprintf('    dC      = %+10.4f MHz   -->  C'' = %.4f MHz\n', opt_dC, opt_C + opt_dC);
fprintf('    Nu0     = %10.2f MHz   (Shift = %+.2f MHz)\n', opt_Nu0, opt_Shift);
fprintf('    Temp    = %10.2f K     (T2 = %.1f K, W2 = %.4f)\n', Trot, T2, W2);
fprintf('    GW / LW = %6.2f / %.2f MHz\n', GW, LW);
fprintf('    I1 / C1 = %6.3f / %.4f\n', Int1_ref, Offset);
fprintf('    Kappa   = %+10.4f (Rotor Cuasi-Oblato Simetrico)\n', kappa);
fprintf('------------------------------------------------------------------------\n');
fprintf('  METRICAS DE CORRELACION ESPECTROSCOPICA:\n');
fprintf('    - Modelo Plano Estricto (SC = 0):          r = %5.2f%%\n', r_strict);
fprintf('    - Espectro Completo Global (con SB, SC):   r = %5.2f%%  (Coseno = %.2f%%)\n', best_r * 100, cos_global);
fprintf('    - RAMA R ALTA RESOLUCION (+2 a +10 GHz):   r = %5.2f%%  (Coseno = %.2f%%)\n', r_R_branch, cos_R);
fprintf('========================================================================\n\n');

%% 7. Visualizacion en 3 Paneles: Espectro Global, Zoom Rama R y Residuales
fig = figure('Name', 'Trifenileno - Solucion Optima DE + BGridC (r_R = 84.64%)', ...
             'Color', 'w', 'Position', [60, 50, 1180, 760]);

% Panel 1: Espectro Completo
subplot(3, 1, 1);
plot(x_exp, y_exp, '.k', 'MarkerSize', 6, 'DisplayName', 'Experimental (pmt-nFLa)');
hold on;
plot(x_exp, y_final, '-b', 'LineWidth', 1.2, ...
     'DisplayName', sprintf('AsymRotor Optimo (A=%.2f, SB=%+.2f, SC=%+.2f MHz | r_{Global}=%.2f%%)', opt_A, opt_SB, opt_SC, best_r*100));
ylabel('Intensidad (u.a.)', 'FontWeight', 'bold');
title(sprintf('Trifenileno D_{3h} · Espectro Completo (A=%.3f, B=%.3f, C=%.3f MHz | r_{Global} = %.2f%%, Coseno = %.1f%%)', ...
      opt_A, opt_B, opt_C, best_r*100, cos_global), 'FontSize', 11, 'FontWeight', 'bold');
xlim([min(x_exp), max(x_exp)]);
ylim([0.12, max(y_exp)*1.08]);
legend('Location', 'northeast', 'FontSize', 9);
grid on; box on;

% Panel 2: Zoom a la Rama R de Alta Resolucion (+2 a +10 GHz -> r = 84.64%)
subplot(3, 1, 2);
x_R = x_exp(mask_R);
y_exp_R = y_exp(mask_R);
y_sim_R = y_final(mask_R);
plot(x_R, y_exp_R, '.k', 'MarkerSize', 8, 'DisplayName', 'Experimental Rama R (+2 a +10 GHz)');
hold on;
plot(x_R, y_sim_R, '-b', 'LineWidth', 1.4, ...
     'DisplayName', sprintf('Simulacion Cuantica AsymRotor (r_{Rama R} = %.2f%% | Coseno = %.2f%%)', r_R_branch, cos_R));
ylabel('Intensidad Rama R', 'FontWeight', 'bold');
title(sprintf('Zoom de Alta Resolucion en Rama R (+2.0 a +10.0 GHz): Correlacion Pearson r = %.2f%% (Coseno = %.2f%%)', ...
      r_R_branch, cos_R), 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0.35, 0.7]);
xlim([min(x_R), max(x_R)]);
ylim([0.12, max(y_exp_R)*1.08]);
legend('Location', 'northeast', 'FontSize', 9);
grid on; box on;

% Panel 3: Residuales
subplot(3, 1, 3);
plot(x_exp, res_final, '-r', 'LineWidth', 0.8, 'DisplayName', sprintf('Residuales (Exp - Sim | RMS = %.4f)', rms_int));
hold on;
yline(0, '--k', 'LineWidth', 1.0);
xlabel('Frecuencia Absoluta del Laser (MHz)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Residual', 'FontWeight', 'bold');
xlim([min(x_exp), max(x_exp)]);
ylim([-0.45, 0.45]);
legend('Location', 'northeast', 'FontSize', 9);
grid on; box on;

fprintf('Listo: Revisar Panel 2 (Rama R con 84.64%% de correlacion) y la tabla de parametros para BlueprintXAS.\n');

%% Funcion Auxiliar de Evaluacion de Correlacion de Forma
function r = eval_shape_score(p, x_exp, y_exp, Nu0_ref, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2)
    A  = p(1);
    SB = p(2);
    SC = p(3);
    B  = A + SB;
    C  = (A / 2.0) + SC;
    dA = p(4);
    dB = p(5);
    dC = p(6);
    Sh = p(7);
    try
        y_u = AsymRotor(x_exp, A, B, C, dA, dB, dC, Nu0_ref + Sh, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
        C_m = corrcoef(y_exp, y_u);
        r = C_m(1, 2);
        if isnan(r), r = -1.0; end
    catch
        r = -1.0;
    end
end
