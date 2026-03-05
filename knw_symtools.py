import numpy as np
import sympy as sp
import symcontroltools as sct

def ptest():
    print("aaa")
    return 0



def addparam(basedict : dict, pdict : dict) -> dict:
    """
    同じキーを持ちタプルを値とする2つの辞書pdict,basedictについて、
    pdictの要素を、basedictの対応するタプルの末尾に追加して返す
    """
    dkeys = basedict.keys()
    for key in dkeys:
        if (str(key) in pdict):
            basedict[key] = sp.Tuple(*basedict[key], pdict.get(str(key)))

def re_aug(mat):
    return mat[:, :-1], mat[:, -1]


def make_sym_canonform(systuple: tuple, mode: str):
    """
    A,B,C,Dのsympy行列を含むタプルを受け取り、標準形を作って返す。
    mode:
    "diag": 対角標準形
    "obs" : 可観測標準形
    delmatnums: すべて1にする行列の番号(主に対角標準形について、B行列かC行列を1に正規化する場合)
    """
    sysdim, iodim = systuple[1].shape
    A = sp.zeros(sysdim, sysdim)
    B = sp.zeros(sysdim, iodim)
    C = sp.zeros(iodim, sysdim)
    D = sp.zeros(iodim, iodim)

    if mode == "diag":
        B = systuple[1]
        C = systuple[2]
        D = systuple[3]
    if mode == "obs":
        B = systuple[1]
        D = systuple[3]

    for i in range(sysdim):
        for j in range(sysdim):
            if mode == "diag":
                if i == j:
                    A[i,j] = systuple[0][i,j]
            if mode == "obs":
                if j == (sysdim - 1):
                    A[i,j] = systuple[0][i,j]
                if (i - j) == 1:
                    A[i,j] = 1
                if (i == 0) and (j == sysdim - 1):
                    C[i,j] = 1
    

    return A,B,C,D



def partial_rref(matrix, limit_col: int):
    """
    指定した列数(limit_col)までで掃き出し法を止める関数
    """
    # 元の行列を書き換えないようコピーを作成
    M = matrix.copy()
    rows, cols = M.shape
    current_row = 0

    # 指定した列数までループ
    for j in range(min(cols, limit_col)):
        if current_row >= rows:
            break

        # 1. ピボット選択
        pivot_row = -1
        for i in range(current_row, rows):
            if M[i, j] != 0:
                pivot_row = i
                break
        
        if pivot_row == -1:
            continue # この列はすべてゼロなのでスキップ

        # 2. 行の入れ替え
        if pivot_row != current_row:
            M.row_swap(current_row, pivot_row)

        # 3. ピボットを1にする（正規化）
        pivot_val = M[current_row, j]
        # 行全体を pivot_val で割る
        M[current_row, :] = M[current_row, :] / pivot_val

        # 4. 掃き出し（他の行のこの列の成分をゼロにする）
        for i in range(rows):
            if i != current_row:
                factor = M[i, j]
                # 行基本変形: row[i] = row[i] - factor * row[current_row]
                if factor != 0:
                    M[i, :] = M[i, :] - factor * M[current_row, :]

        current_row += 1

    return M

def elim_zero_row(M):
    M = sp.Matrix(M).as_mutable()
    rows = M.tolist()
    return sp.Matrix([row for row in rows if any(x != 0 for x in row)])
    #return sp.Matrix([row for row in M.tolist() if not M([row]).is_zero])

def solve_S(aug, vsyms):
    """
    拡大係数行列と解きたい変数のリストを取り、解きたい変数分行標準形にした行列と、解いた変数のdictを返す
    """
    E = aug.echelon_form()
    E = elim_zero_row(E)
    E = partial_rref(E, len(vsyms))
    for i in range(len(vsyms)):
        for j in range(len(vsyms)):
            if E[i,j] == -1:
                E[i,:] = -1*E[i,:]
    ssdict = {}
    if E[:len(vsyms), :len(vsyms)] == sp.eye(len(vsyms)):# s_i = ... の形に解けているか
        for i,v in enumerate(vsyms):
            ssdict[v] = E[i, -1]
    else:
        raise NameError
    return E, ssdict



def systodiag(A, B, C, D):
    egnvals, T = np.linalg.eig(A)
    idx = np.argsort(egnvals)[::-1]
    egnvals = egnvals[idx]
    T = T[:,idx]

    T = T * (np.linalg.inv(T) @ B).reshape(1,T.shape[0])
    Tinv = np.linalg.inv(T)

    Apd = Tinv @ A @ T
    Bpd = Tinv @ B
    Cpd = C @ T
    Dpd = D
    tol = 1e-10
    # 絶対値がtol未満の場所を0にする
    Apd[np.abs(Apd) < tol] = 0.0; Bpd[np.abs(Bpd) < tol] = 0.0; Cpd[np.abs(Cpd) < tol] = 0.0
    return (Apd, Bpd, Cpd, Dpd)

def tf_to_ssmat(G, s):
        G = sp.simplify(sp.expand(G))
        et1 = sct.sym_to_ssobs(G, s)
        return et1.A.T, et1.C.T, et1.B.T, et1.D

