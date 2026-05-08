function ss = zpk_to_observable_canonical(zeros_, poles_, gain)
%ZPK_TO_OBSERVABLE_CANONICAL
% MATLAB port of symcontools.py zpk_to_observable_canonical()
% Pipeline:
%   zpk -> (num,den) -> pad -> controllable canonical -> observable canonical
%
% Returns a struct ss with fields A,B,C,D (complex allowed).



    ss = controllable_to_observable( ...
            zpk_to_controllable_canonical(zeros_, poles_, gain) );
end

% -------------------------------------------------------------------------
% --- make cont and obs ss by zpk (ported) ---
% -------------------------------------------------------------------------

function x = as_1d_complex(x)
% _as_1d_complex(x): None -> [], else -> complex row vector
    if nargin == 0 || isempty(x)
        x = complex([]);  % empty
        return;
    end
    x = complex(x(:).'); % row
end

function [num, den] = zpk_to_tf(zeros_, poles_, gain)
% Convert ZPK to TF polynomials (descending powers).
% Returns (num, den) as complex row vectors.
    z = as_1d_complex(zeros_);
    p = as_1d_complex(poles_);
    k = complex(gain);

    % poly([]) -> 1  (MATLAB behavior)
    num = k * poly(z);   % descending
    den = poly(p);       % descending (monic)
    num = complex(num);
    den = complex(den);
end

function num = pad_num_to_den(num, den)
% Make numerator length equal to denominator length (proper TF expected).
% If deg(num) < deg(den), pad on the LEFT with zeros.
% If equal, keep as-is (D term exists).
    num = complex(num(:).'); % row
    den = complex(den(:).'); % row
    n_den = numel(den);
    if numel(num) < n_den
        num = [zeros(1, n_den - numel(num)), num];
    end
end

function ss = tf_to_controllable_canonical(num, den)
% Build controllable canonical (companion) realization for SISO TF:
%   P(s) = num(s) / den(s), den monic (den(1)=1),
%   num length == den length (pad if needed).
%
% Matches the Python implementation:
%   D = num(1)
%   C = num(2:end) - num(1)*den(2:end)
%   A companion with last row -den(2:end) reversed
%   B = [0;...;0;1]
%
% Returns struct with fields A,B,C,D.

    num = complex(num(:).'); % row
    den = complex(den(:).'); % row

    n = numel(den) - 1;      % order

    % Enforce lengths (minimal handling)
    num = pad_num_to_den(num, den);

    % Companion A (controllable canonical)
    A = complex(zeros(n, n));
    if n >= 2
        A(1:n-1, 2:n) = eye(n-1);
    end
    A(n, :) = -fliplr(den(2:end));

    B = complex(zeros(n, 1));
    B(n, 1) = 1.0;

    D = complex(num(1));  % 1x1 scalar
    C = fliplr(complex(num(2:end) - D * den(2:end))); % 1xn row

    ss = struct('A', A, 'B', B, 'C', C, 'D', D);
end

function ss_o = controllable_to_observable(ss_c)
% Observable canonical form for SISO can be obtained by duality:
%   A_o = A_c^T,  B_o = C_c^T,  C_o = B_c^T,  D_o = D_c
    ss_o = struct();
    ss_o.A = ss_c.A.';
    ss_o.B = ss_c.C.';
    ss_o.C = ss_c.B.';
    ss_o.D = ss_c.D;
end

function ss_c = zpk_to_controllable_canonical(zeros_, poles_, gain)
% Main pipeline:
%   zpk -> (num,den) -> pad -> controllable canonical
    [num, den] = zpk_to_tf(zeros_, poles_, gain);
    num = pad_num_to_den(num, den);

    ss_c = tf_to_controllable_canonical(num, den);

    % Clean up tiny imaginary parts for readability (Python used real_if_close tol=100000)
    tol = 1e5;
    ss_c.A = real_if_close(ss_c.A, tol);
    ss_c.B = real_if_close(ss_c.B, tol);
    ss_c.C = real_if_close(ss_c.C, tol);
    ss_c.D = real_if_close(ss_c.D, tol);
end

function y = real_if_close(x, tol)
% Rough MATLAB equivalent of numpy.real_if_close(x, tol=...).
% If imaginary parts are sufficiently small relative to eps-scale, drop them.
    x = complex(x);

    scale = max(1, max(abs(x(:))));
    thresh = tol * eps(scale);

    if all(abs(imag(x(:))) <= thresh)
        y = real(x);
    else
        y = x;
    end
end