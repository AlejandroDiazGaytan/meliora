function [optParams, stats] = fit_trifenileno_asymrotor(dataFile)
%FIT_TRIFENILENO_ASYMROTOR  Ajuste del espectro LIF de trifenileno con AsymRotor.
%   Evalua la parametrizacion cuasi-oblata con desviaciones simetricas:
%       B = A + SB,   C = A/2 + SC
%   utilizando el modelo:
%       y = I1 * 1000 * AsymRotor(x, A, A+SB, A/2+SC, dA, dB, dC, ...
%                                 Nu0+Shift, Temp, Theta, Phi, GW, LW, ...
%                                 JMax, Unit, SpinAxis, SpinOdd, T2, W2) + C1
%
%   Refinamiento en dos etapas: muestreo diferencial en caja abierta para SC
%   seguido de busqueda en rejilla convergente (3 ciclos de contraccion).

if nargin < 1 || isempty(dataFile)
    dataFile = 'pmt-nFLa';
end

if ~exist(dataFile, 'file')
    error('fit_trifenileno_asymrotor:FileNotFound', ...
        'No se encontro el archivo de datos: %s', dataFile);
end

raw   = load(dataFile);
x_exp = raw(:, 1);
y_exp = raw(:, 2);

% Parametros fijos del perfil instrumental y expansion supersonica
Nu0_ref  = 9.10238e8;     % Origen base en ecuacion BlueprintXAS (910238000 MHz)
Trot     = 3.04;          % Temperatura rotacional primaria (K)
Theta    = 44.8;          % Orientacion dipolar en el plano (grados)
Phi      = 0.0;
GW       = 31.53;         % FWHM Gaussiano (MHz)
LW       = 3.10;          % FWHM Lorentziano (MHz)
JMax     = 50;
Unit     = 1;             % 1 = MHz
SpinAxis = 2;             % Eje b
SpinOdd  = 2.0;
T2       = 12.0;          % Segunda componente termica (K)
W2       = 0.0501;        % Peso de T2
C1       = 0.2211;        % Offset de linea base
I1       = 4.28287;       % Escala global (para I1*1000*AsymRotor)

% 1. Referencia con condicion plana estricta (SB = 0, SC = 0 -> C = A/2)
y_ref0 = AsymRotor(x_exp, 503.0, 503.0, 251.5, 10.0, 8.0, 4.0, ...
                   910237170.0, 6.97, 90.0, 0.0, 37.20, 4.29, 50, 1, 2, 1.0);
R0 = corrcoef(y_exp, y_ref0);
r_planar_strict = R0(1, 2);

% 2. Exploracion de parametros activos: p = [A, SB, SC, dA, dB, dC, Shift]
% Nota: SC se extiende hasta +30 MHz fuera del intervalo [-10, 10] inicial.
p0 = [506.88854, -0.44159, 26.39768, -3.09507, -3.09941, -1.46120, 1915.19];
lb = p0 - [1.5, 0.5, 2.0, 0.5, 0.5, 0.5, 15.0];
ub = p0 + [1.5, 0.5, 2.0, 0.5, 0.5, 0.5, 15.0];

NP   = 12;
nGen = 3;
F    = 0.70;
CR   = 0.85;

rng(42);
pop = repmat(lb, NP, 1) + rand(NP, 7) .* repmat(ub - lb, NP, 1);
pop(1, :) = p0;

scores = zeros(NP, 1);
for i = 1:NP
    scores(i) = shapeScore(pop(i, :), x_exp, y_exp, Nu0_ref, Trot, ...
                           Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
end
[best_r, idx] = max(scores);
best_p = pop(idx, :);

for g = 1:nGen
    for i = 1:NP
        rp = randperm(NP, 3);
        rp(rp == i) = [];
        v = best_p + F * (pop(rp(1), :) - pop(rp(2), :));
        v = max(lb, min(ub, v));
        
        u = pop(i, :);
        jr = randi(7);
        for j = 1:7
            if rand() < CR || j == jr
                u(j) = v(j);
            end
        end
        
        sc_u = shapeScore(u, x_exp, y_exp, Nu0_ref, Trot, ...
                          Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
        if sc_u >= scores(i)
            pop(i, :) = u;
            scores(i) = sc_u;
            if sc_u > best_r
                best_r = sc_u;
                best_p = u;
            end
        end
    end
end

% 3. Contraccion local en rejilla (3 ciclos sobre A, SC, Shift)
dA_grid  = 0.08;
dSC_grid = 0.08;
dSh_grid = 1.50;

for c = 1:3
    vA  = linspace(best_p(1) - dA_grid,  best_p(1) + dA_grid,  3);
    vSC = linspace(best_p(3) - dSC_grid, best_p(3) + dSC_grid, 3);
    vSh = linspace(best_p(7) - dSh_grid, best_p(7) + dSh_grid, 3);
    
    for i1 = 1:3
        for i2 = 1:3
            for i3 = 1:3
                p_try    = best_p;
                p_try(1) = vA(i1);
                p_try(3) = vSC(i2);
                p_try(7) = vSh(i3);
                sc = shapeScore(p_try, x_exp, y_exp, Nu0_ref, Trot, ...
                                Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
                if sc > best_r
                    best_r = sc;
                    best_p = p_try;
                end
            end
        end
    end
    dA_grid  = 0.4 * dA_grid;
    dSC_grid = 0.4 * dSC_grid;
    dSh_grid = 0.4 * dSh_grid;
end

% 4. Evaluacion final del modelo
A_opt  = best_p(1);
SB_opt = best_p(2);
SC_opt = best_p(3);
B_opt  = A_opt + SB_opt;
C_opt  = 0.5 * A_opt + SC_opt;
dA_opt = best_p(4);
dB_opt = best_p(5);
dC_opt = best_p(6);
Sh_opt = best_p(7);
Nu0_opt = Nu0_ref + Sh_opt;

y_unit = AsymRotor(x_exp, A_opt, B_opt, C_opt, dA_opt, dB_opt, dC_opt, ...
                   Nu0_opt, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);

y_calc = (I1 * 1000.0) * y_unit + C1;
resid  = y_exp - y_calc;
rms_y  = sqrt(mean(resid.^2));
kappa  = (2*B_opt - A_opt - C_opt) / (A_opt - C_opt);

% Correlacion por regiones (Global, Rama R completa y ventana de alta S/N)
rel_ghz = (x_exp - Nu0_opt) / 1000.0;
m_R1    = (rel_ghz >= 2.0) & (rel_ghz <= 10.0);
m_R2    = (rel_ghz >= 6.5) & (rel_ghz <= 9.5);

R_r1 = corrcoef(y_exp(m_R1), y_calc(m_R1));
R_r2 = corrcoef(y_exp(m_R2), y_calc(m_R2));
r_R_branch = R_r1(1, 2);
r_R_window = R_r2(1, 2);

% Salida por consola (formato tabla BlueprintXAS)
fprintf('\n----------------------------------------------------------------------\n');
fprintf(' Ajuste AsymRotor - Trifenileno (%s, N = %d puntos)\n', dataFile, numel(x_exp));
fprintf(' Modelo: I1*1000*AsymRotor(x, A, A+SB, A/2+SC, dA, dB, dC, Nu0+Shift, ...)+C1\n');
fprintf('----------------------------------------------------------------------\n');
fprintf(' Parametro      Valor Optimo       Constante Derivada\n');
fprintf('----------------------------------------------------------------------\n');
fprintf(' A              %10.4f MHz     A" = %10.4f MHz\n', A_opt, A_opt);
fprintf(' SB             %+10.4f MHz     B" = %10.4f MHz\n', SB_opt, B_opt);
fprintf(' SC             %+10.4f MHz     C" = %10.4f MHz\n', SC_opt, C_opt);
fprintf(' dA             %+10.4f MHz     A'' = %10.4f MHz\n', dA_opt, A_opt + dA_opt);
fprintf(' dB             %+10.4f MHz     B'' = %10.4f MHz\n', dB_opt, B_opt + dB_opt);
fprintf(' dC             %+10.4f MHz     C'' = %10.4f MHz\n', dC_opt, C_opt + dC_opt);
fprintf(' Shift          %+10.2f MHz     Nu0 = %11.2f MHz\n', Sh_opt, Nu0_opt);
fprintf(' Temp / T2      %6.2f / %4.1f K    W2  = %10.4f\n', Trot, T2, W2);
fprintf(' Theta / Phi    %6.2f / %4.1f deg  Kappa = %+8.4f\n', Theta, Phi, kappa);
fprintf(' GW / LW        %6.2f / %4.2f MHz  SpinOdd = %6.2f (SpinAxis=%d)\n', GW, LW, SpinOdd, SpinAxis);
fprintf(' I1 / C1        %6.3f / %6.4f\n', I1, C1);
fprintf('----------------------------------------------------------------------\n');
fprintf(' Correlacion Pearson (SC = 0, plano estricto) : %.4f (%.2f%%)\n', r_planar_strict, 100*r_planar_strict);
fprintf(' Correlacion Pearson (Global, banda completa) : %.4f (%.2f%%)\n', best_r, 100*best_r);
fprintf(' Correlacion Pearson (Rama R, +2 a +10 GHz)   : %.4f (%.2f%%)\n', r_R_branch, 100*r_R_branch);
fprintf(' Correlacion Pearson (Rama R, +6.5 a +9.5 GHz): %.4f (%.2f%%)\n', r_R_window, 100*r_R_window);
fprintf(' RMS residual de intensidad                   : %.4f\n', rms_y);
fprintf('----------------------------------------------------------------------\n\n');

% Estructuras de salida
optParams = struct('A', A_opt, 'SB', SB_opt, 'SC', SC_opt, ...
                   'B', B_opt, 'C', C_opt, ...
                   'dA', dA_opt, 'dB', dB_opt, 'dC', dC_opt, ...
                   'Nu0', Nu0_opt, 'Shift', Sh_opt, ...
                   'Temp', Trot, 'Theta', Theta, 'Phi', Phi, ...
                   'GW', GW, 'LW', LW, 'I1', I1, 'C1', C1);
stats = struct('r_global', best_r, 'r_R_branch', r_R_branch, ...
               'r_R_window', r_R_window, 'rms', rms_y, 'kappa', kappa);

% Grafica comparativa
figure('Name', 'AsymRotor Fit - Trifenileno', 'Color', 'w', 'Position', [80, 60, 1100, 720]);

subplot(3, 1, 1);
plot(x_exp, y_exp, '.k', 'MarkerSize', 6); hold on;
plot(x_exp, y_calc, '-b', 'LineWidth', 1.1);
ylabel('Intensity (a.u.)');
title(sprintf('Trifenileno (pmt-nFLa) | A=%.3f, B=%.3f, C=%.3f MHz (SB=%+.2f, SC=%+.2f) | r_{global}=%.2f%%', ...
      A_opt, B_opt, C_opt, SB_opt, SC_opt, 100*best_r));
legend('Exp (pmt-nFLa)', 'AsymRotor Fit', 'Location', 'northeast');
xlim([min(x_exp), max(x_exp)]);
ylim([0.12, max(y_exp)*1.08]);
grid on; box on;

subplot(3, 1, 2);
plot(x_exp(m_R1), y_exp(m_R1), '.k', 'MarkerSize', 7); hold on;
plot(x_exp(m_R1), y_calc(m_R1), '-b', 'LineWidth', 1.3);
ylabel('Intensity (R-branch)');
title(sprintf('R-Branch (+2.0 to +10.0 GHz) | r = %.2f%% (Window +6.5 to +9.5 GHz: r = %.2f%%)', ...
      100*r_R_branch, 100*r_R_window));
xlim([min(x_exp(m_R1)), max(x_exp(m_R1))]);
ylim([0.12, max(y_exp(m_R1))*1.08]);
grid on; box on;

subplot(3, 1, 3);
plot(x_exp, resid, '-r', 'LineWidth', 0.8); hold on;
yline(0, '--k');
xlabel('Frequency (MHz)');
ylabel('Residual');
xlim([min(x_exp), max(x_exp)]);
ylim([-0.45, 0.45]);
grid on; box on;

end

function r = shapeScore(p, x_exp, y_exp, Nu0_ref, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2)
A  = p(1);
SB = p(2);
SC = p(3);
B  = A + SB;
C  = 0.5 * A + SC;
dA = p(4);
dB = p(5);
dC = p(6);
Sh = p(7);
try
    y_u = AsymRotor(x_exp, A, B, C, dA, dB, dC, Nu0_ref + Sh, ...
                    Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2);
    R = corrcoef(y_exp, y_u);
    r = R(1, 2);
    if isnan(r), r = -1.0; end
catch
    r = -1.0;
end
end
