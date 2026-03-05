from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, Sequence

import sympy as sp
import os
import numpy as np
import pandas as pd
import control as ct
import io
import json
from pathlib import Path


def dp(ep, capt="", pmode=0, is_expt=False, isptype = False, is_latex = False):
	"""
	sympyで計算した数式をjupyter notebook上で表示する際、latex表現とsrepr()テキストを同時に表示する。
	引数 pmode : 2のときはそのままテキスト出力 
	isptype : 渡された変数の型と(存在すれば)形を出力するか
	depend :
	import sympy as sp
	import os
	"""
	if pmode == 2:
		display(ep)
		return
	pstr = str(capt)# +": "
	if isptype:
		pstr += str(type(ep))
		if hasattr(ep, "shape"):
			pstr += ", " + str(ep.shape)
	#dstr = "$$" + sp.latex(ep) + "$$" + os.linesep + "<details><summary>sympy.srepr()</summary><p>" + sp.srepr(ep) + "</p></details>"
	dstr = "$$" + sp.latex(ep) + "$$" + "<details><summary>sympy.srepr()</summary><p>" + sp.srepr(ep) + "</p></details>"
	if str(capt) != "": 
		display(pstr)
	display(ep)#; print(dstr)
	if is_expt == True:
		print(dstr)
	if is_latex:
		print(sp.latex(ep))
	return

		

def setf(tp, suf, size):
    """
    連番の添え字付きのシンボリック変数を生成するテキストを生成する。
    tpは本体の文字、
    例: {a_{p}}_{0}, ..., {a_{p}}_{2} : setf("a", "p", 3)
    """
    return "{" + tp + "_{"+suf + "}}" + "_{1:"+str(size + 1)+"}"

def get_SISO_sims(dim, suf):
    A = sp.Matrix( a := sp.symbols(setf("a", suf, dim**2)) ).reshape(dim, dim)
    B = sp.Matrix( b := sp.symbols(setf("b", suf, dim)) ).reshape(dim, 1)
    C = sp.Matrix( c := sp.symbols(setf("c", suf, dim)) ).reshape(1, dim)
    D = sp.Matrix( [0] )
    return [(a, b, c), (A, B, C, D)]


def eqs_to_mat(eqList: list, valList: list) -> sp.Matrix:
	"""
	シンボリック式のリストを行列に変換する。
	ex:
	eqList : List
	valList = [ap[1], ap[3], bp[0], bp[1], 1]
	"""
	return sp.Matrix([sp.Matrix([ sp.Poly(eq, valList[:-1]).coeff_monomial(v) for v in valList ]).reshape(1,len(valList)) for eq in list(eqList) ])

def eqs_to_mateqs(eqList: list, valList: list) -> sp.Matrix:
	"""
	aaa
	"""
	mat = eqs_to_mat(eqList, valList)
	matA = mat[:, :-1]
	matX = sp.Matrix(valList[:-1])
	matY = -1*mat[:, -1]
	return ( sp.Eq(sp.MatMul(matA, matX, evaluate=False), matY, evaluate=False), matA, matX, matY )


def getDataDict(csv_order, csv_data: str, isFile: bool):
	if isFile:
		df = pd.read_csv(csv_data, header=None, names=csv_order, skipinitialspace=True, dtype=float)
	else:
		f = io.StringIO(csv_data.strip())
		df = pd.read_csv(f, header=None, names=csv_order, skipinitialspace=True, dtype = float)
	datadict = df.to_dict(orient='list')
	return datadict

def setDataOrder(arg_order, datadict):
	ordered_values = [datadict[key] for key in arg_order]
	val = [ list(d) for d in list(zip(*ordered_values))]
	return val

def get_tf(A, B, C, s):
    tf = C * (s*sp.eye(A.shape[0]) - A).inv() * B
    return tf


def ss_to_ssobs(sys):
	n=sys.A.shape[0]
	T = np.zeros([n,n])
	for i in range(n):
		T[i, (n-1)-i] = 1
	return ct.similarity_transform(sys, T)

def zpk_to_ss_obs(zeros, poles, gain):
	zpksys = ct.zpk(zeros=zeros, poles=poles, gain=gain)
	sssys  = ss_to_ssobs(ct.tf2ss(zpksys))
	return (sssys.A, sssys.B, sssys.C, sssys.D), zpksys

def sym_to_ssobs(tfs, sym):
	nc = np.array(sp.poly(sp.numer(tfs), sym).all_coeffs())
	dc = np.array(sp.poly(sp.denom(tfs), sym).all_coeffs())
	nc = [float(i) for i in nc]
	dc = [float(i) for i in dc]
	edc = ct.tf(*[nc, dc])
	sssys  = ct.tf2ss(edc)
	return ss_to_ssobs(sssys)

def readsacjson(keys):
	script_path = Path(__file__).resolve()
	script_dir = script_path.parent
	json_path = script_dir / "srepr_texts.json"
	with open(json_path, 'r', encoding='utf-8') as f:
		d = json.load(f)
	return [d[k] for k in keys]
	




# --- make cont and obs ss by zpk ---

def _as_1d_complex(x: Iterable[complex] | None) -> np.ndarray:
    if x is None:
        return np.array([], dtype=complex)
    return np.asarray(list(x), dtype=complex).ravel()


def zpk_to_tf(zeros: Sequence[complex] | None,
              poles: Sequence[complex] | None,
              gain: complex = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert ZPK to transfer function polynomials (descending powers).
    Returns (num, den) as complex numpy arrays.
    """
    z = _as_1d_complex(zeros)
    p = _as_1d_complex(poles)
    k = complex(gain)

    # poly([]) -> [1]
    num = k * np.poly(z)  # descending
    den = np.poly(p)      # descending
    return np.asarray(num, dtype=complex), np.asarray(den, dtype=complex)



def pad_num_to_den(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    """
    Make numerator length equal to denominator length (proper TF expected).
    If num degree < den degree, pad on the left with zeros.
    If num degree == den degree, keep as-is (D term exists).
    """

    n_den = den.size
    if num.size < n_den:
        num = np.pad(num, (n_den - num.size, 0), mode="constant")
    return num


def tf_to_controllable_canonical(num: np.ndarray, den: np.ndarray) -> tuple:
    """
    Build controllable canonical (companion) realization for SISO TF:
      P(s) = num(s) / den(s), with den monic (den[0]=1),
      num length == den length (pad if needed).

    This matches the standard tf2ss-style construction:
      D = num[0]
      C = num[1:] - num[0]*den[1:]
      A companion with last row -den[1:][::-1]
      B = [0,0,...,1]^T
    """
    num = np.asarray(num, dtype=complex).ravel()
    den = np.asarray(den, dtype=complex).ravel()

    # Order
    n = den.size - 1

    # Enforce lengths (minimal handling)
    num = pad_num_to_den(num, den)

    # Companion A (controllable canonical)
    A = np.zeros((n, n), dtype=complex)
    A[:-1, 1:] = np.eye(n - 1, dtype=complex)
    A[-1, :] = -den[1:][::-1]

    B = np.zeros((n, 1), dtype=complex)
    B[-1, 0] = 1.0

    D = np.array([[num[0]]], dtype=complex)
    C = (num[1:] - num[0] * den[1:]).reshape(1, n) #厳密にプロパーならnum[0]==0

    return A, B, C, D


def controllable_to_observable(ss: tuple) -> tuple:
    """
    Observable canonical form for SISO can be obtained by duality:
      A_o = A_c^T,  B_o = C_c^T,  C_o = B_c^T,  D_o = D_c
    """
    A, B, C, D = ss[0], ss[1], ss[2], ss[3]
    return A.T, C.T, B.T, D.copy()


def zpk_to_controllable_canonical(zeros: Sequence[complex] | None,
                                poles: Sequence[complex] | None,
                                gain: complex = 1.0) -> tuple:
    """
    Main pipeline:
      zpk -> (num, den) -> normalize -> pad -> controllable canonical
    """
    num, den = zpk_to_tf(zeros, poles, gain)
    num = pad_num_to_den(num, den)

    ss_c = tf_to_controllable_canonical(num, den)

    # Clean up tiny imaginary parts for readability (optional)
    ss_c = (
        np.real_if_close(ss_c[0], tol=100000),
        np.real_if_close(ss_c[1], tol=100000),
        np.real_if_close(ss_c[2], tol=100000),
        np.real_if_close(ss_c[3], tol=100000),
    )
    return ss_c

def zpk_to_observable_canonical(zeros: Sequence[complex] | None,
                                poles: Sequence[complex] | None,
                                gain: complex = 1.0) -> tuple:
    """
    Main pipeline:
      zpk -> controllable canonical -> observable canonical
    """
    return controllable_to_observable(zpk_to_controllable_canonical(zeros=zeros, poles=poles, gain=gain))


def ss_to_mats(ss):
	return(ss.A, ss.B, ss.C, ss.D)

# ------------------------------
# Example usage
# ------------------------------
if __name__ == "__main__":
    # Example: P(s) = 2 (s+1) / ((s+2)(s+3)) = 2*(s - (-1))/((s - (-2))(s - (-3)))
    zeros = [-1.0]
    poles = [-2.0, -3.0]
    k = 2.0

    ss = zpk_to_observable_canonical(zeros, poles, k)
    print("A=\n", ss.A)
    print("B=\n", ss.B)
    print("C=\n", ss.C)
    print("D=\n", ss.D)







