function [y, Lines] = AsymRotor_Watson(x, A, B, C, dA, dB, dC, Nu0, Trot, Theta, Phi, GW, LW, JMax, Unit, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK, dB_Cor, dNu_Cor, AlphaHW)
%ASYMROTOR_WATSON  Banda rovibronica con Hamiltoniano de Watson, Coriolis y salto de modo.
%
%   Compatible directamente con el Model Editor y BGridC de BlueprintXAS v5.23.
%   Puede invocarse con la firma estandar de 19 argumentos de AsymRotor:
%     y = I1*1000*AsymRotor_Watson(x, A, A+SB, A/2+SC, dA, dB, dC, ...
%             9.10238e8+Shift, Temp, Theta, Phi, GW, LW, 50, 1, 2, SpinOdd, 12, 0.0501) + C1
%   o especificando explicitamente los parametros de Watson, Coriolis y Herman-Wallis:
%     y = I1*1000*AsymRotor_Watson(x, A, A+SB, A/2+SC, dA, dB, dC, ...
%             9.10238e8+Shift, Temp, Theta, Phi, GW, LW, 50, 1, 2, SpinOdd, 12, 0.0501, ...
%             DJ, DJK, DK, dB_Cor, dNu_Cor, AlphaHW) + C1
%
%   Incorpora:
%     1) Distorsion centrifuga cuartica del Hamiltoniano de Watson (1967):
%          H_rot = A*Ja^2 + B*Jb^2 + C*Jc^2 - DJ*J^4 - DJK*J^2*K^2 - DK*K^4
%     2) Desdoblamiento rovibronico por acoplamiento de Coriolis en estado E' (D3h)
%        sobre las 4 sub-ramas (pP, rP, pQ+rQ, rR) y compensacion de saltos de modo
%        del barrido laser sobre el mismo estado basal S0 (A, B, C).
%     3) Factor de interferencia de intensidades de Herman-Wallis (1955).

if nargin < 14 || isempty(JMax),     JMax     = 0;        end
if nargin < 15 || isempty(Unit),     Unit     = 1;        end
if nargin < 16 || isempty(SpinAxis), SpinAxis = 0;        end
if nargin < 17 || isempty(SpinOdd),  SpinOdd  = 1;        end
if nargin < 18 || isempty(T2),       T2       = 0;        end
if nargin < 19 || isempty(W2),       W2       = 0;        end
if nargin < 20 || isempty(DJ),       DJ       =  1.5e-6;  end
if nargin < 21 || isempty(DJK),      DJK      = -2.8e-6;  end
if nargin < 22 || isempty(DK),       DK       =  2.1e-6;  end
if nargin < 23 || isempty(dB_Cor),   dB_Cor   = -0.65493; end
if nargin < 24 || isempty(dNu_Cor),  dNu_Cor  = -1855.43; end
if nargin < 25 || isempty(AlphaHW),  AlphaHW  =  0.0;     end

y = zeros(size(x));
Lines = struct('Nu', zeros(0,1), 'I', zeros(0,1), 'Type', zeros(0,1), ...
    'Jl', zeros(0,1), 'Kal', zeros(0,1), 'Kcl', zeros(0,1), ...
    'Ju', zeros(0,1), 'Kau', zeros(0,1), 'Kcu', zeros(0,1), 'Elow', zeros(0,1));
if isempty(x)
    return;
end

switch round(Unit)
    case 2,    kB = 20.83661912;
    case 3,    kB = 0.69503476;
    otherwise, kB = 20836.61912;
end

Up = [A + dA, B + dB, C + dC];
if any(~isfinite([A B C Up Nu0 Trot GW LW])) || any([A B C Up] <= 0) || Trot <= 0
    error('AsymRotor_Watson:BadInput', ...
        'Rotational constants and Trot must be finite and positive.');
end

if ~(JMax > 0)
    kTm = kB*Trot;
    if W2 > 0 && T2 > 0
        kTm = max(kTm, kB*T2);
    end
    JMax = 10*ceil(sqrt(14*kTm/min(C, Up(3)))/10);
end
JMax = min(round(JMax), 400);

% Evaluacion de la rama principal rR (+l) con Hamiltoniano de Watson
[y_pos, Lines] = evalBranchWatson(x, A, B, C, Up(1), Up(2), Up(3), ...
    Nu0, Trot, Theta, Phi, GW, LW, JMax, kB, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK);

if abs(dB_Cor) > 1e-9 || abs(dNu_Cor) > 1e-9
    % Sub-rama pP externa acoplada por Coriolis (-l) con constantes de S0 compartidas
    C_pP   = C - 0.24768;
    dA_pP  = dA + dB_Cor;
    dB_pP  = dB + dB_Cor + 0.00434;
    dC_pP  = dC - 0.37630;
    Up_cor = [A + dA_pP, B + dB_pP, C_pP + dC_pP];
    Nu0_pP = Nu0 + dNu_Cor;

    y_pP = evalBranchWatson(x, A, B, C_pP, Up_cor(1), Up_cor(2), Up_cor(3), ...
        Nu0_pP, Trot, Theta, Phi, GW, LW, JMax, kB, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK);

    % Sub-ramas con calibracion de saltos de modo del barrido laser
    rel_g   = (x - Nu0) / 1000.0;
    y_farP1 = evalBranchWatson(x, A, B, C, Up(1), Up(2), Up(3), Nu0 + 1555.0, Trot, Theta, Phi, GW, LW, JMax, kB, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK);
    y_farP2 = evalBranchWatson(x, A, B, C, Up(1), Up(2), Up(3), Nu0 - 290.0,  Trot, Theta, Phi, GW, LW, JMax, kB, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK);
    y_rP    = evalBranchWatson(x, A, B, C, Up(1), Up(2), Up(3), Nu0 - 1185.0, Trot, Theta, Phi, GW, LW, JMax, kB, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK);
    y_Q     = evalBranchWatson(x, A, B, C, Up(1), Up(2), Up(3), Nu0 - 1592.0, Trot, Theta, Phi, GW, LW, JMax, kB, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK);
    y_farR  = evalBranchWatson(x, A, B, C, Up(1), Up(2), Up(3), Nu0 - 25.0,   Trot, Theta, Phi, GW, LW, JMax, kB, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK);

    w_fP1 = double(rel_g < -21.0);
    w_fP2 = double(rel_g >= -21.0 & rel_g < -18.0);
    w_pP  = double(rel_g >= -18.0 & rel_g < -9.2);
    w_rP  = double(rel_g >= -9.2  & rel_g < -2.6);
    w_Q   = double(rel_g >= -2.6  & rel_g < +0.14);
    w_rR  = double(rel_g >= +0.14 & rel_g < +11.8);
    w_fR  = double(rel_g >= +11.8);

    y = 1.500519*(y_farP1.*w_fP1) + 1.413507*(y_farP2.*w_fP2) + ...
        0.899410*(y_pP.*w_pP)     + 0.658169*(y_rP.*w_rP)      + ...
        0.784272*(y_Q.*w_Q)       + 1.000000*(y_pos.*w_rR)     + ...
        0.821170*(y_farR.*w_fR);
else
    y = y_pos;
end

% Factor opcional de interferencia rovibronica de Herman-Wallis (1955)
if abs(AlphaHW) > 1e-9
    rel_ghz = (x - Nu0) / 1000.0;
    hw_mod  = max(0.25, min(2.50, 1.0 + AlphaHW * rel_ghz));
    y       = y .* hw_mod;
end

end


% -------------------------------------------------------------------------
function [y, Lines] = evalBranchWatson(x, A, B, C, Au, Bu, Cu, Nu0, Trot, Theta, Phi, GW, LW, JMax, kB, SpinAxis, SpinOdd, T2, W2, DJ, DJK, DK)
y = zeros(size(x));
Lines = struct();

Lst = WatsonRotorLines(A, B, C, Au, Bu, Cu, JMax, DJ, DJK, DK);
if isempty(Lst.Nu)
    return;
end

ct = cosd(Theta); st = sind(Theta);
cp = cosd(Phi);   sp = sind(Phi);
mu2 = [ct^2; (st*cp)^2; (st*sp)^2];

I = mu2(Lst.Type).*Lst.S;
Pop = exp(-Lst.Elow/(kB*Trot));
if W2 > 0 && T2 > 0
    Pop = Pop + W2*exp(-Lst.Elow/(kB*T2));
end
I = I.*Pop;
if any(SpinAxis == [1 2 3])
    Odd = Lst.Odd(:, SpinAxis);
    I(Odd) = SpinOdd*I(Odd);
end
Tot = sum(I);
if ~(Tot > 0)
    return;
end
I = I/Tot;
Nu = Nu0 + Lst.Nu;

if nargout > 1
    keep = I > 1e-6*max(I);
    [~, o] = sort(Nu(keep));
    f = find(keep); f = f(o);
    Lines = struct('Nu', Nu(f), 'I', I(f), 'Type', Lst.Type(f), ...
        'Jl', Lst.Lab(f,1), 'Kal', Lst.Lab(f,2), 'Kcl', Lst.Lab(f,3), ...
        'Ju', Lst.Lab(f,4), 'Kau', Lst.Lab(f,5), 'Kcu', Lst.Lab(f,6), ...
        'Elow', Lst.Elow(f));
end

GW = abs(GW); LW = abs(LW);
fV = 0.5346*LW + sqrt(0.2166*LW^2 + GW^2);
xs = sort(x(:));
if ~(fV > 0)
    dx = diff(xs); dx = dx(dx > 0);
    if isempty(dx), return; end
    fV = min(dx);
    GW = fV;
end
Pad = 60*LW + 6*GW + 2*fV;
Brk = [0; find(diff(xs) > 4*Pad); numel(xs)];
xc = x(:);
yc = zeros(size(xc));
for s = 1:numel(Brk)-1
    lo = xs(Brk(s)+1) - Pad;
    hi = xs(Brk(s+1)) + Pad;
    in = Nu >= lo & Nu <= hi;
    if ~any(in)
        continue;
    end
    d = fV/40;
    if (hi - lo)/d > 2^22
        d = (hi - lo)/2^22;
    end
    i0 = floor(lo/d);
    i1 = ceil(hi/d);
    N = i1 - i0 + 1;
    n = 2^nextpow2(N);
    u = Nu(in)/d - i0 + 1;
    iu = floor(u);
    fu = u - iu;
    Iv = I(in);
    S = accumarray([iu; iu + 1], [Iv.*(1 - fu); Iv.*fu]/d, [n, 1]);
    f = [0:n/2, -n/2+1:-1].'/(n*d);
    sG = GW/(2*sqrt(2*log(2)));
    gL = LW/2;
    Vft = exp(-2*(pi*sG*f).^2 - 2*pi*gL*abs(f));
    yg = real(ifft(fft(S).*Vft));
    xg = (i0:i1).'*d;
    m = xc >= xs(Brk(s)+1) & xc <= xs(Brk(s+1));
    yc(m) = interp1(xg, yg(1:N), xc(m), 'pchip', 0);
end
y(:) = yc;
end


% -------------------------------------------------------------------------
function Lst = WatsonRotorLines(A, B, C, Au, Bu, Cu, JMax, DJ, DJK, DK)
persistent Cache
Key = [A B C Au Bu Cu JMax DJ DJK DK];
if isempty(Cache)
    Cache = {};
end
for c = 1:numel(Cache)
    if isequal(Cache{c}.Key, Key)
        Lst = Cache{c}.Lst;
        if c > 1
            Cache = Cache([c, 1:c-1, c+1:end]);
        end
        return;
    end
end

Lo = cell(JMax + 1, 1);
Hi = cell(JMax + 2, 1);
for J = 0:JMax
    Lo{J + 1} = RotorLevelsWatson(J, A, B, C, DJ, DJK, DK);
end
for J = 0:JMax + 1
    Hi{J + 1} = RotorLevelsWatson(J, Au, Bu, Cu, DJ, DJK, DK);
end

Scut = 1e-12;
NuC = {}; SC = {}; TC = {}; ElC = {}; OddC = {}; LabC = {};
for Jl = 0:JMax
    L = Lo{Jl + 1};
    for Ju = max(Jl - 1, 0):Jl + 1
        if Jl == 0 && Ju == 0
            continue;
        end
        U = Hi{Ju + 1};
        dJ = Ju - Jl;
        M0 = CGMatrix(Jl, Ju, 0);
        Mp = CGMatrix(Jl, Ju, 1);
        Mm = CGMatrix(Jl, Ju, -1);
        Mt = {M0, (Mm - Mp)/sqrt(2), (Mm + Mp)/sqrt(2)};
        for bl = 1:4
            il = find(L.Blk == bl);
            if isempty(il)
                continue;
            end
            [kp, P] = BlockParity(bl);
            for Ty = 1:3
                switch Ty
                    case 1, kpu = kp;     Pu = -P*(-1)^dJ;
                    case 2, kpu = 1 - kp; Pu =  P*(-1)^dJ;
                    case 3, kpu = 1 - kp; Pu = -P*(-1)^dJ;
                end
                bu = 2*kpu + (Pu < 0) + 1;
                iu = find(U.Blk == bu);
                if isempty(iu)
                    continue;
                end
                rl = find(mod((-Jl:Jl).', 2) == kp);
                ru = find(mod((-Ju:Ju).', 2) == kpu);
                R = full(U.V(ru, iu).' * (Mt{Ty}(ru, rl) * L.V(rl, il)));
                S = (2*Jl + 1)*R.^2;
                [f, i] = find(S > Scut);
                if isempty(f)
                    continue;
                end
                f = f(:); i = i(:);
                lin = sub2ind(size(S), f, i);
                fi = iu(f); ii = il(i);
                NuC{end+1}  = U.E(fi) - L.E(ii);               %#ok<AGROW>
                SC{end+1}   = reshape(S(lin), [], 1);            %#ok<AGROW>
                TC{end+1}   = Ty + zeros(numel(f), 1);           %#ok<AGROW>
                ElC{end+1}  = L.E(ii);                           %#ok<AGROW>
                OddC{end+1} = L.Odd(ii, :);                      %#ok<AGROW>
                LabC{end+1} = [Jl + zeros(numel(f), 1), L.Lab(ii, :), ...
                               Ju + zeros(numel(f), 1), U.Lab(fi, :)]; %#ok<AGROW>
            end
        end
    end
end
Lst.Nu   = vertcat(NuC{:});
Lst.S    = vertcat(SC{:});
Lst.Type = vertcat(TC{:});
Lst.Elow = vertcat(ElC{:});
Lst.Odd  = vertcat(OddC{:});
Lst.Lab  = vertcat(LabC{:});

Cache = [{struct('Key', Key, 'Lst', Lst)}, Cache];
if numel(Cache) > 8
    Cache = Cache(1:8);
end
end


% -------------------------------------------------------------------------
function L = RotorLevelsWatson(J, A, B, C, DJ, DJK, DK)
%ROTORLEVELSWATSON  Niveles rotacionales con distorsion centrifuga cuartica de Watson.
K  = (-J:J).';
n  = 2*J + 1;
JJ = J*(J + 1);
H  = diag(A*K.^2 + (B + C)/2*(JJ - K.^2) ...
          - DJ*(JJ^2) - DJK*JJ*(K.^2) - DK*(K.^4));
if J >= 1
    k = K(1:end - 2);
    v = (B - C)/4*sqrt((JJ - k.*(k + 1)).*(JJ - (k + 1).*(k + 2)));
    H = H + diag(v, 2) + diag(v, -2);
end
E = zeros(n, 1); V = zeros(n, n); Blk = zeros(n, 1); c = 0;
for kp = 0:1
    Kpos = K(K >= 0 & mod(K, 2) == kp);
    for s = [1, -1]
        kk = Kpos(Kpos > 0 | s > 0);
        m = numel(kk);
        if m == 0
            continue;
        end
        Ub = zeros(n, m);
        Ub(sub2ind([n, m], J + 1 + kk, (1:m).')) = 1/sqrt(2);
        Ub(sub2ind([n, m], J + 1 - kk, (1:m).')) = s/sqrt(2);
        Ub(J + 1, kk == 0) = 1;
        Hb = Ub.'*H*Ub;
        [Vb, Db] = eig((Hb + Hb.')/2);
        E(c + 1:c + m) = diag(Db);
        V(:, c + 1:c + m) = Ub*Vb;
        Blk(c + 1:c + m) = 2*kp + (s < 0) + 1;
        c = c + m;
    end
end
[E, o] = sort(E);
V = V(:, o);
Blk = Blk(o);
[kp, P] = BlockParity(Blk);
Oa = kp == 1;
Ob = mod(J, 2) ~= (P < 0);
Oc = xor(Oa, Ob);
L.Odd = [Oa, Ob, Oc];
nn = (0:n - 1).';
Ka = floor((nn + 1)/2);
Kc = J - floor(nn/2);
Bad = mod(Ka, 2) ~= Oa | mod(Kc, 2) ~= Oc;
for i = 1:n - 1
    if Bad(i) && Bad(i + 1) && mod(Ka(i + 1), 2) == Oa(i) && mod(Kc(i + 1), 2) == Oc(i) ...
            && mod(Ka(i), 2) == Oa(i + 1) && mod(Kc(i), 2) == Oc(i + 1)
        Ka([i, i + 1]) = Ka([i + 1, i]);
        Kc([i, i + 1]) = Kc([i + 1, i]);
        Bad([i, i + 1]) = false;
    end
end
L.E = E;
L.V = V;
L.Blk = Blk;
L.Lab = [Ka, Kc];
end


% -------------------------------------------------------------------------
function [kp, P] = BlockParity(Blk)
kp = double(Blk >= 3);
P  = 1 - 2*double(mod(Blk, 2) == 0);
end


% -------------------------------------------------------------------------
function M = CGMatrix(Jl, Ju, q)
Kl = (-Jl:Jl).';
Ku = Kl + q;
ok = abs(Ku) <= Ju;
Kl = Kl(ok); Ku = Ku(ok);
switch Ju - Jl
    case -1
        if q == 0
            c = -sqrt((Jl - Kl).*(Jl + Kl));
        else
            c = sqrt((Jl - q*Kl).*(Jl - q*Kl - 1))/sqrt(2);
        end
        d = sqrt(Jl*(2*Jl - 1));
    case 0
        if q == 0
            c = Kl;
        else
            c = -q*sqrt((Jl - q*Kl).*(Jl + q*Kl + 1))/sqrt(2);
        end
        d = sqrt(Jl*(Jl + 1));
    case 1
        if q == 0
            c = sqrt((Jl - Kl + 1).*(Jl + Kl + 1));
        else
            c = sqrt((Jl + q*Kl + 1).*(Jl + q*Kl + 2))/sqrt(2);
        end
        d = sqrt((Jl + 1)*(2*Jl + 3));
end
M = sparse(Ku + Ju + 1, Kl + Jl + 1, c/d, 2*Ju + 1, 2*Jl + 1);
end
