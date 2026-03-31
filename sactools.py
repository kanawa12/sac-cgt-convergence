import numpy as np
import sympy as sp
from sympy.physics.control import StateSpace
import symcontools as sct

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



def get_syms_tuple(m):
    return sorted(m.free_symbols, key=sp.default_sort_key)

def build_ss_plant_and_refmodel(pdim: int, mdim:int, pmode: str, mmode: str):
    _, sp_bases = sct.get_SISO_sims(pdim,"p")
    ss_p = StateSpace(*sp_bases)
    _, sm_bases = sct.get_SISO_sims(mdim,"m")
    ss_m = StateSpace(*sm_bases)
    app = sct.make_sym_canonform(ss_p, pmode, True)
    amp = sct.make_sym_canonform(ss_m, mmode, True)
    return ss_p, ss_m, app, amp, {"Ap": get_syms_tuple(app.A), "Bp": get_syms_tuple(app.B), "Am": get_syms_tuple(amp.A), "Bm": get_syms_tuple(amp.B)}



def make_smats(pdim: int, rdim: int) -> sp.BlockMatrix:
    sasyms = sct.makesyms("s", "x", rdim, rdim)
    sbsyms = sct.makesyms("s", "u", pdim, 1)
    Kxsyms = sct.makesyms("k", "x", 1, rdim)
    Kusyms = sct.makesyms("k", "u", 1, 1)
    return sp.BlockMatrix([[sasyms, sbsyms], [Kxsyms, Kusyms]]), {"Sx": get_syms_tuple(sasyms), "Su": get_syms_tuple(sbsyms), "Kx": get_syms_tuple(Kxsyms), "Ku": get_syms_tuple(Kusyms)}



def sym_setting(pdim, rdim, maindict):
    maindict["sp_bases"], maindict["sm_bases"], maindict["sp_obs"], maindict["sm_obs"] = build_ss_plant_and_refmodel(pdim, pdim, "obs", "obs")
    maindict["sskk_bases"] = make_smats(pdim, rdim)
    return maindict

def make_CGTEqs(sblockmat: sp.BlockMatrix, sp_obs: StateSpace, sm_obs: StateSpace):
    """
    CGTのS行列(blockmatrix)とプラント、規範モデルのsympy statespaceを受け取り、CGT条件式を計算して返す
    この時点では式は整理しない
    """
    Aps = sct.ss_to_blockmat(sp_obs)
    Ams = sct.ss_to_blockmat(sm_obs)
    sx = sblockmat.blocks[0,0]; n = sx.rows
    scoeff = sp.BlockMatrix([
    [sx, sp.ZeroMatrix(n, 1)],
    [sp.ZeroMatrix(1, n), sp.Identity(1)]])
    expr = Aps*sblockmat - scoeff*Ams
    KAS = sp.Tuple(Aps, sblockmat, scoeff, Ams)
    return expr.as_explicit(), KAS

def solve_CGTEq(cgteq, vals: tuple) -> dict:
    """
    cgt方程式の各成分のリストと変数のタプルをとり、解いた結果の辞書を返す
    """
    ssmat, ssmaty = sp.linear_eq_to_matrix(cgteq, vals)
    E, ssdict = solve_S(sp.Matrix([[ssmat, ssmaty]]), vals)
    return ssdict


def solve_CGTEq2(expr, ssdict, apsyms: tuple, bpsyms: tuple, kxsyms: tuple, kusyms: tuple):
    """
    CGT条件式とSx,Suを解いた辞書を受け取り、プラントパラメータやCGT解(Kx,Ku)について整理した式を返す
    """
    exprs2 = elim_zero_row(list(sp.simplify(expr.subs(ssdict)))) # Sx,Suを消去し整理したCGT条件式
    absyms: tuple = apsyms + bpsyms
    abmat = sp.Matrix(absyms).reshape(len(absyms), 1) #未知数ベクトル
    kxkusyms: tuple = kxsyms + kusyms
    kmat  = sp.Matrix(kxkusyms).reshape(len(kxkusyms), 1)

    mats = sp.linear_eq_to_matrix(exprs2, absyms)
    th = partial_rref(sp.Matrix([[*mats]]), len(apsyms) )
    #th[-1, :] /= th[-1, -1] # 最後の行を正規化
    abans = sp.Tuple(th[:, :-1], abmat, th[:, -1])
    # 方程式を求める
    reeq = sp.expand(th[:, :-1]@abmat - th[:, -1])
    keqmats = sp.linear_eq_to_matrix(reeq, kxkusyms)
    kans = sp.Tuple(keqmats[0], kmat, keqmats[1])
    return reeq, abans, kans

def CGT_def_and_solve(pdim, rdim):
    _, _, ss_plant, ss_refmodel, symdict = build_ss_plant_and_refmodel(pdim, rdim, "obs", "obs")
    # 制御対象と規範モデルのsym StateSpaceの作成 symdictはシンボリック変数のリスト

    SMat, sdict = make_smats(pdim, rdim)
    # SAC,CGTで出てくるS_x, S_u, K_x, K_uの定義

    symdict = symdict | sdict

    expr, kas = make_CGTEqs(SMat, ss_plant, ss_refmodel)
    # CGT方程式の作成

    ssdict = solve_CGTEq(expr, svals:=symdict["Sx"] + symdict["Su"])
    # CGT方程式をS行列について解き、解を辞書の形で返す

    eq, abans, kans = solve_CGTEq2(expr, ssdict, symdict["Ap"], symdict["Bp"], symdict["Kx"], symdict["Ku"])
    # CGT方程式にS行列を代入して消去し、(Ap,Bp), (Kx,Ku) それぞれについてまとめる
    resultdict = {
        "plant": ss_plant, "refmodel": ss_refmodel, "SMatrix": SMat, "CGT_expr": expr, "KAS": kas,
        "ssdict": ssdict, "plainEq":eq, "apbpEq": abans, "kxkuEq": kans
    }
    return resultdict
