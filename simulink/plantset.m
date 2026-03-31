clear
% DEMO: zpk のリスト(zeros, poles, gain)から
% 1) 伝達関数係数を作る
% 2) 可制御標準形(Controllable Canonical Form; CCF)を作る
% 3) CCF を使って可観測標準形(Observable Canonical Form; OCF)へ変換
%
% ※この OCF は「CCF の転置」として作るので、A の符号規約は
%   D(s)=s^n + a_{n-1}s^{n-1}+...+a0 のとき A の最終列に -[a0;...;a_{n-1}] が入ります。
%   ご提示の a_{in} が正符号で欲しい場合は、最後に A(:,end) の符号を反転してください。

Sp1 = zpk2struct([-1, -3], [2, -5+1j, -5-1j], 3) %[output:35aec887]

Sm1 = zpk2struct([-1, -3]    , [-1.5, -7+3j, -7-3j], 1) %[output:47b05b2d]
Sm2 = zpk2struct([-1, -3]    , [-2, -4, -4.5], 2) %[output:712b66b8]
Sm3 = zpk2struct([-1.1, -3]  , [-2, -4, -4.5], 2) %[output:98e2b69b]
Sm4 = zpk2struct([-1.1, -3.2], [-2, -4, -4.5], 2) %[output:22fd56f8]

% clipboard('copy', str)
Sp = Sp1; Sm = Sm4; S = Sm;
refname = "{S_{m}}_{4}";
reffullname = "Ref.Model 4";

%%


out = sim("SAC_noise_ident_simu_3.slx"); %[output:8545790f]

texs = ssblock2latex(S)
tez  = zpklist2latex(reffullname, S)
few = mat4x3_to_latex(out.cgt.signals.values, refname)
tgb = est_to_latex(out.est.signals.values, refname)

clipboard('copy', few);
clipboard('copy', tgb);



%%

zs = 0.0:0.1:1.0; zs
Css = []
for i = zs
    disp(i);
    Sm = zpk2struct([-i, -3] , [-1.5, -7+3j, -7-3j], 1);
    out = sim("SAC_noise_ident_simu_3.slx");
    Css = [Css, out.cgt.signals.values(:, end)];
end

cs = Css;
disp(Css)

%save('zerodif.mat','zs','cs', '-v7')
%%
Sm = zpk2struct([-0.7, -3]    , [-1.5, -7+3j, -7-3j], 1)
%%


A = rand(4,3);
B = rand(4,2);

A_last = A(:, end);      % 4x1（最終列）
C = [B, A_last];         % hstack（横に連結）
% つまり C は 4x3
%%
g = out.glogs;

%u_t = g.err.Time;
t = g.err.Time;
u   = g.u_m.Data;
yp = g.y_p.Data;
%k_t = g.kxkulog.Time;
kxku   = g.kxkulog.Data;   % 2501x4
ke = g.k_e.Data;
err = g.err.Data;

save('simlog2.mat','t','u','kxku','ke',"err","yp",'-v7')

%%


function latex = est_to_latex(M, suf)

    assert(all(size(M) == [3 3]), 'Input must be 3x3 matrix');

    % 行ラベル（固定）
    labels = {
        '{a_p}_{11;' + suf + '}'
        '{a_p}_{12;' + suf + '}'
        '{a_p}_{13;' + suf + '}'
    };

    latex = strings(4,1);

    for i = 1:3
        c1 = fmt(M(i,1), 5);   % 高精度
        c2 = fmt(M(i,2), 5);   % 高精度
        c3 = fmt(M(i,3), 3);   % 粗め

        latex(i) = sprintf('%s & %s & %s & %s \\\\', ...
                           labels{i}, c1, c2, c3);
    end

    latex(end) = latex(end) + " \hline";
    latex = strjoin(latex, newline);
end

function latex = mat4x3_to_latex(M, suf)

    assert(all(size(M) == [4 3]), 'Input must be 4x3 matrix');

    % 行ラベル（固定）
    labels = {
        '{k_x}_{11;' + suf + '}'
        '{k_x}_{12;' + suf + '}'
        '{k_x}_{13;' + suf + '}'
        '{k_u}_{11;' + suf + '}'
    };

    latex = strings(4,1);

    for i = 1:4
        c1 = fmt(M(i,1), 5);   % 高精度
        c2 = fmt(M(i,2), 5);   % 高精度
        c3 = fmt(M(i,3), 3);   % 粗め

        latex(i) = sprintf('%s & %s & %s & %s \\\\', ...
                           labels{i}, c1, c2, c3);
    end

    latex(end) = latex(end) + " \hline";
    latex = strjoin(latex, newline);
end


function s = fmt(x, sig)

    if x == 0
        s = '0';
        return
    end

    e = floor(log10(abs(x)));
    m = x / 10^e;

    m = round(m, sig - 1);

    if abs(e) <= 2
        s = sprintf(['%.' num2str(sig-1) 'g'], x);
    else
        s = sprintf(['%.' num2str(sig-1) 'g \\times 10^{%d}'], m, e);
    end
end






function s_s = zpk2struct(z, p, k)
s_s.z = z; s_s.p = p; s_s.k = k;
ss = zpk_to_observable_canonical(z, s_s.p, s_s.k);
s_s.A = ss.A; s_s.B = ss.B; s_s.C = ss.C; s_s.D = ss.D;
end


function tex = ssblock2latex(s_s)
%SSBLOCK2LATEX  Return LaTeX string of a block state-space matrix [A|B; C|D]
A = s_s.A; B = s_s.B; C = s_s.C; D = s_s.D;
fmt = "%.6g";  % change if you want (e.g., "%.4g", "%.8f")

n = size(A,1);
m = size(B,2);
p = size(C,1);

% alignment like "cc|c"
align = [repmat('c',1,n) '|' repmat('c',1,m)];

num2s = @(x) string(sprintf(fmt, x));

% build body rows
rows = strings(n+p,1);

Mtop = [A B];
for i = 1:n
    parts = strings(1,n+m);
    for j = 1:n+m
        parts(j) = num2s(Mtop(i,j));
    end
    rows(i) = strjoin(parts, " & ");
end

Mbot = [C D];
for i = 1:p
    parts = strings(1,n+m);
    for j = 1:n+m
        parts(j) = num2s(Mbot(i,j));
    end
    rows(n+i) = strjoin(parts, " & ");
end

bodyTop = strjoin(rows(1:n), " \\ " + newline);
bodyBot = strjoin(rows(n+1:end), " \\ " + newline);



tex = "$$" + newline + ...
" S_{p}:\left[\begin{array}{c|c}" + newline + ...
" A_{p} & B_{p} \\ \hline" + newline + ...
" C_{p} & D_{p}" + newline + ...
" \end{array}\right]=" + newline + ...
" \left[\begin{array}{" + string(align) + "}" + newline + ...
bodyTop + " \\ \hline" + newline + ...
bodyBot + newline + ...
" \end{array}\right]" + newline + ...
"$$";

tex = char(tex);  % return as char (single string), like typical LaTeX snippet
end


function tex = zpklist2latex(sysNames, s_s)
%ZPKLIST2LATEX  Build LaTeX table for multiple SISO ZPK sets
%
% Inputs (assumed valid; no extra error handling):
%   sysNames  : string/cellstr, N systems (e.g., ["Plant","Ref. Model 1",...])
%   zerosList : cell array (1xN or Nx1), each element is numeric vector of zeros
%   polesList : cell array (1xN or Nx1), each element is numeric vector of poles
%   gains     : numeric vector length N
%
% Output:
%   tex : char (single LaTeX string)


zerosList = s_s.z; polesList = s_s.p; gains = s_s.k;

fmt = "%.6g";   % number format

N = numel(gains);

rows = strings(N,1);
for k = 1:N
    zstr = vec2latex(zerosList(k,:), fmt);
    pstr = vec2latex(polesList(k,:), fmt);
    gstr = string(sprintf(fmt, gains(k)));

    rows(k) = "\text{" + string(sysNames(k)) + "} & " + zstr + " & " + pstr + " & " + gstr + " \\";
end

body = strjoin(rows, newline);


tex0 = "$$" + newline + ...
"\begin{array}{c|ccc} \hline" + newline + ...
"\text{System} & \text{Zeros} & \text{Poles} & \text{Gain} \\ \hline" + newline + ...
body + newline + ...
"\hline" + newline + ...
"\end{array}" + newline + ...
"$$";


tex = body + " \hline" + newline;

tex = char(tex);
end


function s = vec2latex(v, fmt)
% Convert numeric vector to LaTeX-friendly comma-separated string.
% Special-case: show complex conjugate pair as "a\pm bi" if detected.

v = v(:).'; % row

if isempty(v)
    s = "";
    return
end

% Detect one conjugate pair [a+bi, a-bi] (any order), show a\pm bi
if numel(v) == 2 && ~isreal(v(1)) && ~isreal(v(2))
    if abs(real(v(1)) - real(v(2))) < 1e-12 && abs(imag(v(1)) + imag(v(2))) < 1e-12
        a = real(v(1));
        b = abs(imag(v(1)));
        s = real2s(a, fmt) + "\pm " + real2s(b, fmt) + "i";
        return
    end
end

parts = strings(1,numel(v));
for i = 1:numel(v)
    parts(i) = cplx2s(v(i), fmt);
end
s = strjoin(parts, ", ");
end


function s = real2s(x, fmt)
s = string(sprintf(fmt, x));
end

function s = cplx2s(x, fmt)
if isreal(x)
    s = real2s(x, fmt);
else
    a = real(x);
    b = imag(x);
    if b >= 0
        s = real2s(a, fmt) + "+" + real2s(b, fmt) + "i";
    else
        s = real2s(a, fmt) + real2s(b, fmt) + "i"; % b already has "-"
    end
end
end

%%





% --- 実行 ---
[num, den] = zpk_list_to_tf_poly(Sp_z, Sp_p, Sp_k);
[Aco,Bco,Cco,Dco] = tf_poly_to_ccf(num, den);
[A,B,C,D]         = ccf_to_ocf(Aco,Bco,Cco,Dco);

disp("=== TF coefficients ===");
disp("num = "); disp(num);
disp("den = "); disp(den);

disp("=== Controllable canonical (Aco,Bco,Cco,Dco) ===");
disp(Aco); disp(Bco); disp(Cco); disp(Dco);

disp("=== Observable canonical (A,B,C,D) ===");
disp(A); disp(B); disp(C); disp(D);



function [num, den] = zpk_list_to_tf_poly(z, p, k)
% z, p: ベクトル（複素数可）
% k: スカラー
% num, den: s の降冪順係数（MATLAB の poly と同じ順）
%
% G(s) = k * prod(s - z_i) / prod(s - p_i)

den = poly(p);            % monic
num = k * poly(z);

% 係数を実数に寄せたい場合（数値誤差で虚部が極小のとき）
num = real_if_close(num);
den = real_if_close(den);
end


function [A,B,C,D] = tf_poly_to_ccf(num, den)
% 伝達関数 N(s)/D(s) から可制御標準形(Controllable companion form)を構成
% den: s^n + a_{n-1}s^{n-1}+...+a0（monic に正規化）
%
% 状態方程式は
% x1dot = x2
% ...
% x_{n-1}dot = x_n
% x_ndot = -a0 x1 - a1 x2 - ... - a_{n-1} x_n + u
%
% B = [0 ... 0 1]^T（入力は最後の状態に入る）
%
% num が den と同次数の場合は直達項 D を含める（多項式割り算で分離）
% ※不適切(分子次数 > 分母次数)はエラー

% row vector 化
num = num(:).';
den = den(:).';

% 正規化（den の先頭係数を 1 に）
if den(1) == 0
    error("Denominator leading coefficient is zero.");
end
num = num / den(1);
den = den / den(1);

n = length(den) - 1;   % 次数

% 不適切（proper でない）チェック
deg_num = poly_degree(num);
deg_den = poly_degree(den);
if deg_num > deg_den
    error("Improper transfer function: deg(num) > deg(den). This code assumes proper TF.");
end

% 分子/分母の多項式割り算で直達項 D と厳密にプロパーな残りを分離
[q, r] = deconv(num, den);   % num = q*den + r

% q は次数差ぶんの多項式。proper なので q は高々定数のはず
q = trim_leading_zeros(q);
if isempty(q)
    D = 0;
elseif length(q) == 1
    D = q(1);
else
    % ここに来たら improper に近い（数値誤差も含む）
    error("Unexpected polynomial part (q has degree > 0).");
end

% r の次数は < n。状態空間の C を作るために長さ n へパディング（降冪順）
r = trim_leading_zeros(r);
if isempty(r)
    r = 0;
end
r_desc = pad_poly_desc(r, n);          % [b_{n-1} ... b0]（降冪順, 長さ n）
b_asc  = fliplr(r_desc);               % [b0 ... b_{n-1}]（昇冪順, 長さ n）

% den = [1 a_{n-1} ... a0]（降冪順）
a_desc = den(2:end);                   % [a_{n-1} ... a0]
a_asc  = fliplr(a_desc);               % [a0 ... a_{n-1}]

% --- 可制御標準形（superdiagonal に 1, 最終行に -a_asc） ---
A = zeros(n,n);
A(1:n-1, 2:n) = eye(n-1);
A(n, :)       = -a_asc;

B = zeros(n,1);
B(n) = 1;

% C は [b0 ... b_{n-1}] - D*[a0 ... a_{n-1}]
C = (b_asc - D * a_asc);

end


function [Ao, Bo, Co, Do] = ccf_to_ocf(Ac, Bc, Cc, Dc)
% 「CCF を用いて OCF を作る」：転置系の関係を使う
%
% (Ac,Bc,Cc,Dc) が可制御標準形なら、
% (Ao,Bo,Co,Do) = (Ac', Cc', Bc', Dc) は可観測標準形になります。
%
% この Ao は
% - subdiagonal に 1（ご提示の形）
% - 最終列が（分母係数由来で） -[a0; a1; ...; a_{n-1}]
%   という形になります。

Ao = Ac.';
Bo = Cc.';
Co = Bc.';
Do = Dc;
end


% ===== ユーティリティ =====

function d = poly_degree(p)
p = trim_leading_zeros(p);
if isempty(p)
    d = -inf;
else
    d = length(p) - 1;
end
end

function p = trim_leading_zeros(p)
% 先頭（高次側）の 0 を落とす（完全に 0 なら空に）
p = p(:).';
idx = find(abs(p) > 0, 1, 'first');
if isempty(idx)
    p = [];
else
    p = p(idx:end);
end
end

function p_desc = pad_poly_desc(p_desc, n)
% p_desc: 降冪順係数（長さ <= n）を、長さ n に左側（高次側）ゼロ詰め
p_desc = p_desc(:).';
if length(p_desc) > n
    % 本来ここは起きない（r の次数は < n のはず）
    p_desc = p_desc(end-n+1:end);
end
p_desc = [zeros(1, n - length(p_desc)), p_desc];
end

function x = real_if_close(x)
% 虚部が極小なら実数化
if ~isreal(x)
    if max(abs(imag(x))) < 1e-12 * max(1, max(abs(real(x))))
        x = real(x);
    end
end
end




%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":36.5}
%---
%[output:35aec887]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"Sp1","value":"    z: [-1 -3]\n    p: [2.0000 + 0.0000i -5.0000 + 1.0000i -5.0000 - 1.0000i]\n    k: 3\n    A: [3×3 double]\n    B: [3×1 double]\n    C: [0 0 1]\n    D: 0\n"}}
%---
%[output:47b05b2d]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"Sm1","value":"    z: [-1 -3]\n    p: [-1.5000 + 0.0000i -7.0000 + 3.0000i -7.0000 - 3.0000i]\n    k: 1\n    A: [3×3 double]\n    B: [3×1 double]\n    C: [0 0 1]\n    D: 0\n"}}
%---
%[output:712b66b8]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"Sm2","value":"    z: [-1 -3]\n    p: [-2 -4 -4.5000]\n    k: 2\n    A: [3×3 double]\n    B: [3×1 double]\n    C: [0 0 1]\n    D: 0\n"}}
%---
%[output:98e2b69b]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"Sm3","value":"    z: [-1.1000 -3]\n    p: [-2 -4 -4.5000]\n    k: 2\n    A: [3×3 double]\n    B: [3×1 double]\n    C: [0 0 1]\n    D: 0\n"}}
%---
%[output:22fd56f8]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"Sm4","value":"    z: [-1.1000 -3.2000]\n    p: [-2 -4 -4.5000]\n    k: 2\n    A: [3×3 double]\n    B: [3×1 double]\n    C: [0 0 1]\n    D: 0\n"}}
%---
%[output:8545790f]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"プログラムの割り込み (Ctrl-C) が検出されました。"}}
%---
