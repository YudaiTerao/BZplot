#!/usr/bin/env python

import re
import sys
from pathlib import Path
from typing import Optional
import numpy as np
from scipy.constants import *
from BZplot.bzplot import cube_path

def bohrang(bohr: float) -> float:
    #bohr: 5.29*10^-11m
    return bohr * physical_constants["Bohr radius"][0]*(10**10)

def _ibravError():
    print("ibravcell: lattice_constants not found")
    sys.exit(0)

def ibravcell(
    ibrav: np.int64,
    A: np.float64 = -1,
    B: np.float64 = -1,
    C: np.float64 = -1,
    cosAB: np.float64 = -2,
    cosAC: np.float64 = -2,
    cosBC: np.float64 = -2
):
    """
    ibravに対応した格子ベクトルを出力する

    Args:
        ibrav(int): ibravの番号
        A, B, C(float): 格子定数
        cosAB, cosAC, cosBC(float): 格子ベクトル間の角度
    """
    # cubic P (sc)
    if ibrav == 1:
        if A <= 0 or \
        not (B == C == -1 and cosAB == cosAC == cosBC == -2): _ibravError()
        cell = np.array([[ 1,  0,  0],
                         [ 0,  1,  0],
                         [ 0,  0,  1]]) * A
    # cubic F (fcc)
    elif ibrav == 2:
        if A <= 0 or \
        not (B == C == -1 and cosAB == cosAC == cosBC == -2): _ibravError()
        cell = np.array([[-1,  0,  1],
                         [ 0,  1,  1],
                         [-1,  1,  0]]) * A / 2
    # cubic I (bcc)
    elif ibrav == 3: 
        if A <= 0 or \
        not (B == C == -1 and cosAB == cosAC == cosBC == -2): _ibravError()
        cell = np.array([[ 1,  1,  1],
                         [-1,  1,  1],
                         [-1, -1,  1]]) * A / 2
    # cubic I (bcc), more symmetric axis
    elif ibrav == -3:
        if A <= 0 or \
        not (B == C == -1 and cosAB == cosAC == cosBC == -2): _ibravError()
        cell = np.array([[-1,  1,  1],
                         [ 1, -1,  1],
                         [ 1,  1, -1]]) * A / 2
    # Hexagonal and Trigonal P
    elif ibrav == 4:
        if A <= 0 or C <= 0 or \
        not (B == -1 and cosAB == cosAC == cosBC == -2): _ibravError()
        cell = np.array([[   1,            0,   0],
                         [-1/2, np.sqrt(3)/2,   0],
                         [   0,            0, C/A]]) * A
    # Trigonal R, 3fold axis c
    elif ibrav == 5:
        if A <= 0 or not B == -1: _ibravError()
        if 0 <= cosAB < 1 and ( cosAC == -2 or cosAC == cosAB) \
                         and ( cosBC == -2 or cosBC == cosAB):cosG = cosAB
        elif 0 <= cosAC < 1 and ( cosAB == -2 or cosAB == cosAC) \
                           and ( cosBC == -2 or cosBC == cosAC):cosG = cosAC
        elif 0 <= cosBC < 1 and ( cosAB == -2 or cosAB == cosBC) \
                           and ( cosAC == -2 or cosAC == cosBC):cosG = cosBC
        elif 0 <= C < 1 :cosG = C
        else: _ibravError()
        tx = np.sqrt((1 - cosG)/2)
        ty = np.sqrt((1 - cosG)/6)
        tz = np.sqrt((1 + 2 * cosG)/3)
        cell = np.array([[  tx,    -ty,  tz],
                         [   0, 2 * ty,  tz],
                         [ -tx,    -ty,  tz]]) * A
    # Trigonal R, 3fold axis <111>
    elif ibrav == -5:
        if A <= 0 or not B == -1: _ibravError()
        if 0 <= cosAB < 1 and ( cosAC == -2 or cosAC == cosAB) \
                         and ( cosBC == -2 or cosBC == cosAB):cosG = cosAB
        elif 0 <= cosAC < 1 and ( cosAB == -2 or cosAB == cosAC) \
                           and ( cosBC == -2 or cosBC == cosAC):cosG = cosAC
        elif 0 <= cosBC < 1 and ( cosAB == -2 or cosAB == cosBC) \
                           and ( cosAC == -2 or cosAC == cosBC):cosG = cosBC
        elif 0 <= C < 1 :cosG = C
        else: _ibravError()
        tx = np.sqrt((1 - cosG)/2)
        ty = np.sqrt((1 - cosG)/6)
        tz = np.sqrt((1 + 2 * cosG)/3)
        A_dash = A / np.sqrt(3)
        u = tz - 2 * np.sqrt(2) * ty
        v = tz + np.sqrt(2) * ty
        cell = np.array([[ u, v, v],
                         [ v, u, v],
                         [ v, v, u]]) * A_dash
    # Tetragonal P (st)
    elif ibrav == 6:
        if A <= 0 or C <= 0 or \
        not (B == -1 and cosAB == cosAC == cosBC == -2): _ibravError()
        cell = np.array([[   1,    0,    0],
                         [   0,    1,    0],
                         [   0,    0,  C/A]]) * A
    # Tetragonal I (bst)
    elif ibrav == 7:
        if A <= 0 or C <= 0 or \
        not (B == -1 and cosAB == cosAC == cosBC == -2): _ibravError()
        cell = np.array([[   1,   -1,  C/A],
                         [   1,    1,  C/A],
                         [  -1,   -1,  C/A]]) * A / 2
    # Orthorhombic P
    elif ibrav == 8:
        if A <= 0 or B <= 0 or C <= 0 or \
        not cosAB == cosAC == cosBC == -2: _ibravError()
        cell = np.array([[ A,  0,  0],
                         [ 0,  B,  0],
                         [ 0,  0,  C]])
    # Orthorhombic base-centered(bco)
    elif ibrav == 9:
        if A <= 0 or B <= 0 or C <= 0 or \
        not cosAB == cosAC == cosBC == -2: _ibravError()
        cell = np.array([[ A/2,  B/2,    0],
                         [-A/2,  B/2,    0],
                         [   0,    0,    C]])
    # Orthorhombic base-centered(bco)( as 9 alternate description )
    elif ibrav == -9:
        if A <= 0 or B <= 0 or C <= 0 or \
        not cosAB == cosAC == cosBC == -2: _ibravError()
        cell = np.array([[ A/2, -B/2,    0],
                         [ A/2,  B/2,    0],
                         [   0,    0,    C]])
    # Orthorhombic one-face base-centered A-type
    elif ibrav == 91:
        if A <= 0 or B <= 0 or C <= 0 or \
        not cosAB == cosAC == cosBC == -2: _ibravError()
        cell = np.array([[   A,    0,    0],
                         [   0,  B/2, -C/2],
                         [   0,  B/2,  C/2]])
    # Orthorhombic face-centered
    elif ibrav == 10:
        if A <= 0 or B <= 0 or C <= 0 or \
        not cosAB == cosAC == cosBC == -2: _ibravError()
        cell = np.array([[ A/2,    0,  C/2],
                         [ A/2,  B/2,    0],
                         [   0,  B/2,  C/2]])
    # Orthorhombic face-centered
    elif ibrav == 11:
        if A <= 0 or B <= 0 or C <= 0 or \
        not cosAB == cosAC == cosBC == -2: _ibravError()
        cell = np.array([[ A/2,  B/2,  C/2],
                         [-A/2,  B/2,  C/2],
                         [-A/2, -B/2,  C/2]])
    # Monoclinic P, unique axis c
    elif ibrav == 12:
        if A <= 0 or B <= 0 or C <= 0 or not -1 < cosAB < 1 or \
        not cosAC == cosBC == -2: _ibravError()
        sinAB = np.sqrt( 1 - cosAB**2 )
        cell = np.array([[         A,         0, 0],
                         [ B * cosAB, B * sinAB, 0],
                         [         0,         0, C]])
    # Monoclinic P, unique axis b
    elif ibrav == -12:
        if A <= 0 or B <= 0 or C <= 0 or not -1 < cosAC < 1 or \
        not cosAB == cosBC == -2: _ibravError()
        sinAC = np.sqrt( 1 - cosAC**2 )
        cell = np.array([[         A,   0,         0],
                         [         0,   B,         0],
                         [ C * cosAC,   0, C * sinAC]])
    # Monoclinic base-centered, unique axis c
    elif ibrav == 13:
        print("ibrav = 13 not support")
#        if A <= 0 or B <= 0 or C <= 0 or not -1 < cosAB < 1: _ibravError()
#        sinAB = np.sqrt( 1 - cosAB**2 )
#        cell = np.array([[       A/2,         0, -C/2],
#                         [ B * cosAB, B * sinAB,    0],
#                         [       A/2,         0,  C/2]])
    # Monoclinic base-centered, unique axis b
    elif ibrav == -13:
        print("ibrav = -13 not support")
#        if A <= 0 or B <= 0 or C <= 0 or not -1 < cosAC < 1: _ibravError()
#        sinAC = np.sqrt( 1 - cosAC**2 )
#        cell = np.array([[       A/2, B/2,         0],
#                         [      -A/2, B/2,         0],
#                         [ C * cosAC,   0, C * sinAC]])
    # Triclinic
    elif ibrav == 14:
        if A <= 0 or B <= 0 or C <= 0 \
        or not -1 < cosAB < 1 \
        or not -1 < cosAC < 1 \
        or not -1 < cosBC < 1: _ibravError()
        sinAB = np.sqrt( 1 - cosAB**2 )
        cell = np.array([[ A, 0, 0],
                         [ B * cosAB, B * sinAB, 0],
                         [ C * cosAC, \
                           C * (cosBC - cosAC * cosAB) / sinAB, \
                           C * np.sqrt(1 + 2 * cosBC * cosAC * cosAB \
                           - cosBC**2 - cosAC**2 - cosAB**2) / sinAB]])
    else: 
        print("ibrav is wrong or not support number")
        sys.exit(0)
    return cell

class BZ_input:
    """
    inputfileを読み込み、オブジェクト変数として格納する

    Attributes:
        cell(array): 実空間の格子ベクトル(3, 3)
        kcell(array): 波数空間の逆格子ベクトル(3, 3)
        kpath(array): kpathの座標(, 3)
        kpath_name(List): kpathのラベル
        atom_frac(array): 原子座標(, 3)
        atom_name(List): 原子のラベル
    """

    def __init__(self, filepath: Optional[Path] = None):
        """
        Args:
            filepath(Path): inputfile path. default is None.
        """
        self.cell = np.zeros((3, 3), dtype = np.float64)
        self.kcell = np.zeros((3, 3), dtype = np.float64)
        self.kpath, self.atom_frac = np.empty(0), np.empty(0)
        self.kpath_name, self.atom_name = [], []

        if filepath == None:
            self.cell, self.kcell = self._test_config()
        elif filepath.exists() == False:
            print("BZ_input: inputfile not found")
            sys.exit(0)
        elif re.match(r'^.*scf[\._]in$', str(filepath)) or \
             re.match(r'^.*nscf[\._]in$', str(filepath)):
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self._read_scf_in(filepath)
        elif re.match(r'^.*[\._]win$', str(filepath)):
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self._read_win(filepath)
        elif re.match(r'^.*wt[\._]in$', str(filepath)) or \
             re.match(r'^.*weyl[\._]in$', str(filepath))or \
             re.match(r'^.*bz[\._]in$', str(filepath)):
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self._read_wt_in(filepath)
        else : 
            print("BZ_input: filename is wrong")
            sys.exit(0)

    def _test_config(self):
        """
        filepathがNoneのとき。動作確認用
        """
        cell = ibravcell(3, A=1.0)
        kcell = 2*pi*(np.linalg.inv(cell).T)
        return cell, kcell

    def _read_scf_in(self, nscf_in_path: Path):
        """
        scf.inやnscf.inを読み込む

        Args:
            scf_in_path(Path): nscf.in or scf.in file path.
        """
        ibrav = 200
        A, B, C, cosAB, cosAC, cosBC = -1, -1, -1, -2, -2, -2
        celldm1, celldm2, celldm3, celldm4, celldm5, celldm6 = -2, -2, -2, -2, -2, -2
        kpath, atom_frac = np.empty(0), np.empty(0)
        kpath_name, atom_name = [], []
        BZ_base, kp_method = "", ""

        with open(nscf_in_path, 'r') as fn:
            #改行,行頭行末スペース,行末のカンマ,コロン,タブ,コメント行を消去
            lines_with_comment = [ re.sub('^!.*$', '', re.sub('^\s+', '', re.sub('[\,\.]?\s+$', '', l))) for l in fn.readlines() ]
            lines = [ l.split('!')[0] for l in lines_with_comment ]

        for i, line in enumerate(lines):
            # lwはlineを小文字に統一し、コメント,空白を削除したもの
            lw = re.sub('\s', '', line).lower().split('!')[0]
            # ibravは必ず整数, 整数でなければ読み込まない
            if re.match(r'^ibrav=-?[0-9]+$', lw):
                ibrav=int(line.split("=")[1])
            #a, b, c, celldmは必ず実数, 実数でなければ読み込まない
            elif re.match(r'^a=-?[0-9]+\.?[0-9]*$', lw):
                A = float(line.split("=")[1])
            elif re.match(r'^b=-?[0-9]+\.?[0-9]*$', lw):
                B = float(line.split("=")[1])
            elif re.match(r'^c=-?[0-9]+\.?[0-9]*$', lw):
                C = float(line.split("=")[1])
            elif re.match(r'^cosab=-?[0-9]+\.?[0-9]*$', lw):
                cosAB = float(line.split("=")[1])
            elif re.match(r'^cosac=-?[0-9]+\.?[0-9]*$', lw):
                cosAC = float(line.split("=")[1])
            elif re.match(r'^cosbc=-?[0-9]+\.?[0-9]*$', lw):
                cosBC = float(line.split("=")[1])
            elif re.match(r'^celldm\(1\)=-?[0-9]+\.?[0-9]*$', lw):
                celldm1 = float(line.split("=")[1])
            elif re.match(r'^celldm\(2\)=-?[0-9]+\.?[0-9]*$', lw):
                celldm2 = float(line.split("=")[1])
            elif re.match(r'^celldm\(3\)=-?[0-9]+\.?[0-9]*$', lw):
                celldm3 = float(line.split("=")[1])
            elif re.match(r'^celldm\(4\)=-?[0-9]+\.?[0-9]*$', lw):
                celldm4 = float(line.split("=")[1])
            elif re.match(r'^celldm\(5\)=-?[0-9]+\.?[0-9]*$', lw):
                celldm5 = float(line.split("=")[1])
            elif re.match(r'^celldm\(6\)=-?[0-9]+\.?[0-9]*$', lw):
                celldm6 = float(line.split("=")[1])

            elif re.match(r'^cell_parameters.+$', lw):
                BZ_base = re.sub(r'^cell_parameters', '', lw)
                cell = np.array([ l.split()[:3] for l in lines[i+1:i+4] ], dtype=np.float64)
            elif re.match(r'^nat=[0-9]+$', lw):
                atom_num = int(line.split("=")[1])
            elif re.match(r'^atomic_positions.+$', lw):
                atom_base = re.sub(r'^atomic_positions', '', lw)
                atom_name = [ l.split()[0] for l in lines[i+1:i+1+atom_num]]
                atom_frac = np.array([ l.split()[1:4] for l in lines[i+1:i+1+atom_num]], dtype=np.float64)
            elif re.match(r'^k_points.+$', lw):
                kp_method = re.sub(r'^k_points', '', lw)
                if 'tpiba_b' in kp_method or 'crystal_b' in kp_method:
                    kp_num = int(lines[i+1])
                    kpath = np.array([ l.split()[:3] for l in lines[i+2:i+2+kp_num] ], dtype=np.float64)
                    for j in range(i+2, i+2+kp_num):
                        if len(lines_[j].split('!')) == 2:
                            kpath_name.append(lines_with_comment[j].split('!')[1])

        #----- cellの単位変更 -----#
        #cell : 実空間の基底ベクトル,単位はangstromに統一
        #cellに関する部分は必ず出力に必要なため、errorの条件分岐を書いている
        if celldm1 != -2: A = bohrang(celldm1)
        if celldm2 != -2: 
            assert celldm1 != -2, "celldm(1) not found"
            B = A * celldm2
        if celldm3 != -2: 
            assert celldm1 != -2, "celldm(1) not found"
            C = A * celldm3

        if ibrav in [0]:
            if 'alat' in BZ_base: 
                if A <= 0 or not (B == C == -1 \
                   and cosAB == cosAC == cosBC == -2) : _ibravError()
                else: cell = cell * A
            elif 'bohr'     in BZ_base: cell = bohrang(cell)
            elif 'angstrom' in BZ_base: pass
            else: 
                print("read_scf_in: CELL_PARAMETERS is wrong")
                sys.exit(0)
        # 角度指定なし
        elif ibrav in [1, 2, 3, -3, 4, 6, 7, 8, 9, -9, 91, 10, 11]:
            cell = ibravcell(ibrav, A = A, B = B, C = C, \
                             cosAB = cosAB, cosAC = cosAC, cosBC = cosBC)
        # 角度指定あり
        # cosAB, cosAC, cosBCを使った角度指定
        elif ibrav in [5, -5, 12, -12, 14] and \
            (celldm4 == celldm5 == celldm6 == -2):
            cell = ibravcell(ibrav, A = A, B = B, C = C, \
                             cosAB = cosAB, cosAC = cosAC, cosBC = cosBC)
        # celldm()を使った角度指定
        elif ibrav in [5, -5, 12, -12, 14] and \
            (cosAB == cosAC == cosBC == -2):
            if ibrav in [5, -5, 12] and (celldm5 == celldm6 == -2):
                cell = ibravcell(ibrav, A = A, B = B, C = C, cosAB = celldm4)
            elif ibrav in [-12] and (celldm4 == celldm6 == -2):
                cell = ibravcell(ibrav, A = A, B = B, C = C, cosAC = celldm5)
            elif ibrav in [14]:
                cell = ibravcell(ibrav, A = A, B = B, C = C, \
                            cosAB = celldm6, cosAC = celldm5, cosBC = celldm4)
            else: _ibravError()
        elif ibrav in [13, -13]:
            #ibrav=13は計算が面倒なのでまだ非対応
            cell = ibravcell(ibrav)
        else:
            print("read_scf_in: ibrav not found or not support number")
            sys.exit(0)

        #----- kcellの計算 -----#
        #kcell : K空間の基底ベクトル
        #kcellは実空間の基底ベクトルの逆行列の転置*2pi
        kcell = 2 * pi * (np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        if 'tpiba_b' in kp_method:
            kpath = kpath * 2 * pi / A
        elif 'crystal_b' in kp_method:
            kpath = np.matmul(kpath, kcell)

        #----- atom_fracの単位変換 -----#
        if   'alat' in atom_base:
            atom_frac = np.matmul( atom_frac * A, np.linalg.inv(cell) )
        elif 'bohr' in atom_base:
            atom_frac = np.matmul( bohrang(atom_frac), np.linalg.inv(cell) )
        elif 'angstrom' in atom_base:
            atom_frac = np.matmul( atom_frac, np.linalg.inv(cell) )
        elif 'crystal' in atom_base: pass

        return cell, kcell, kpath, kpath_name, atom_frac, atom_name

    def _read_win(self, win_path: Path):
        """
        .win fileを読み込む

        Args:
            win_path(Path): win file path.
        """
        kpath_name, atom_name = [], []
        kpath, atom_frac = np.empty(0), np.empty(0)
        BZ_base, atom_base = "", ""
        cell_flag = 0

        with open(win_path, 'r') as fw:
            #改行,行頭行末スペース,タブ,コメントを消去
            lines = [ re.sub('^!.*$', '', re.sub('^\s+', '', re.sub('\s+$', '', l))).split('!')[0] for l in fw.readlines() ]

        for i, line in enumerate(lines):
            # lwはlineを小文字に統一したもの
            lw = line.lower()

            if re.match(r'^begin\s+unit_cell_cart$', lw):
                BZ_base = lines[i+1]
                cell = np.array([ l.split() for l in lines[i+2:i+5] ], \
                                dtype=np.float64)
                if cell.shape == (3, 3): cell_flag = 1
            elif re.match(r'^begin\s+kpoint_path$', lw):
                j, kp = 1, []
                while re.match(r'^end\s+kpoint_path$', lines[i+j].lower()) is None:
                    kpath_name.append(lines[i+j].split()[0])
                    kp.append(lines[i+j].split()[1:4])
                    j = j + 1
                    if i + j >= len(lines): 
                        print("read_win: \'end kpoint_path\' not found")
                        sys.exit(0)
                kpath_name.append(lines[i+j-1].split()[-4])
                kp.append(lines[i+j-1].split()[-3:])
                kpath = np.array(kp, dtype=np.float64)
            elif re.match(r'^begin\s+atoms_frac$', lw) or \
                 re.match(r'^begin\s+atoms_cart$', lw):
                atom_base = lw.split('_')[1]
                j, af = 1, []
                while re.match(r'^end\s+atoms_[a-z]{4}$', lines[i+j].lower()) is None:
                    atom_name.append(lines[i+j].split()[0])
                    af.append(lines[i+j].split()[1:4])
                    j = j + 1
                    if i + j >= len(lines): 
                        print("read_win: \'end atoms_\' not found")
                        sys.exit(0)
                atom_frac = np.array(af, dtype=np.float64)

        #----- cellの単位変更 -----#
        if   'ang'  in BZ_base: pass
        elif 'bohr' in BZ_base: cell = bohrang(cell)
        if ( 'ang' not in BZ_base and 'bohr' not in BZ_base ) or cell_flag == 0:
            print("read_win: unit_cell_cart is wrong")
            sys.exit(0)

        #----- kcellの計算 -----#
        kcell = 2*pi*(np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        if len(kpath) != 0 : kpath = np.matmul(kpath, kcell)

        #----- atom_fracの単位変更 -----#
        if   'frac' in atom_base : pass
        elif 'cart' in atom_base :
            atom_frac = np.matmul( atom_frac, np.linalg.inv(cell) )

        return cell, kcell, kpath, kpath_name, atom_frac, atom_name

    def _read_wt_in(self, wt_in_path: Path):
        """
        WannierToolsのinputのwt.inやweyl.pyのinputのweyl.inを読み込む
        BZplot用のinputのbz.inもこの形式

        Args:
            win_path(Path): win file path.
        """
        ibrav = 200
        A, B, C, cosAB, cosAC, cosBC = -1, -1, -1, -2, -2, -2
        celldm4, celldm5, celldm6 = -2, -2, -2
        kpath_name, atom_name = [], []
        kpath, atom_frac = np.empty(0), np.empty(0)

        with open(wt_in_path, 'r') as fw:
            #改行,行頭行末スペース,タブ,コメントを消去
            lines = [ re.sub('^!.*$', '', re.sub('^\s+', '', re.sub('\s+$', '', l))).split('!')[0] for l in fw.readlines() ]

        for i, line in enumerate(lines):
            # lwはlineを小文字に統一したもの
            lw = line.lower()

            if re.match(r'^lattice$', lw):
                BZ_base = lines[i+1].lower()
                if 'ibrav' in BZ_base:
                    j = 1
                    while i + j < len(lines):
                        lw = lines[i+j].lower()
                        if re.match(r'^ibrav=-?[0-9]+$', lw):
                            ibrav=int(line.split("=")[1])
                        #a, b, c, celldmは必ず実数, 実数でなければ読み込まない
                        elif re.match(r'^a=-?[0-9]+\.?[0-9]*$', lw):
                            A = float(line.split("=")[1])
                        elif re.match(r'^b=-?[0-9]+\.?[0-9]*$', lw):
                            B = float(line.split("=")[1])
                        elif re.match(r'^c=-?[0-9]+\.?[0-9]*$', lw):
                            C = float(line.split("=")[1])
                        elif re.match(r'^cosab=-?[0-9]+\.?[0-9]*$', lw):
                            cosAB = float(line.split("=")[1])
                        elif re.match(r'^cosac=-?[0-9]+\.?[0-9]*$', lw):
                            cosAC = float(line.split("=")[1])
                        elif re.match(r'^cosbc=-?[0-9]+\.?[0-9]*$', lw):
                            cosBC = float(line.split("=")[1])
                        elif re.match(r'^celldm\(1\)=-?[0-9]+\.?[0-9]*$', lw):
                            A = bohrang(float(line.split("=")[1]))
                        elif re.match(r'^celldm\(2\)=-?[0-9]+\.?[0-9]*$', lw):
                            B = bohrang(float(line.split("=")[1]))
                        elif re.match(r'^celldm\(3\)=-?[0-9]+\.?[0-9]*$', lw):
                            C = bohrang(float(line.split("=")[1]))
                        elif re.match(r'^celldm\(4\)=-?[0-9]+\.?[0-9]*$', lw):
                            celldm4 = float(line.split("=")[1])
                        elif re.match(r'^celldm\(5\)=-?[0-9]+\.?[0-9]*$', lw):
                            celldm5 = float(line.split("=")[1])
                        elif re.match(r'^celldm\(6\)=-?[0-9]+\.?[0-9]*$', lw):
                            celldm6 = float(line.split("=")[1])
                        else: break
                        j = j + 1
                else: cell = np.array([ l.split() for l in lines[i+2:i+5] ], \
                                        dtype=np.float64)
            if re.match(r'^atom_positions$', lw):
                atom_num = int(lines[i+1])
                atom_base = lines[i+2].lower()
                atom_name = [ l.split()[0] for l in lines[i+3:i+3+atom_num] ]
                atom_frac = np.array([ l.split()[1:4] for l in lines[i+3:i+3+atom_num] ], dtype=np.float64)

            if re.match(r'^kpath_bulk$', lw):
                kp_num = int(lines[i+1])
                kpath_name = [ l.split()[0] for l in lines[i+2:i+2+kp_num] ]
                kpath_name.append(line[i+1+kp_num].split()[4])
                kp = [ l.split()[1:4] for l in lines[i+2:i+2+kp_num] ]
                kp.append(line[i+1+kp_num].split()[5:8])
                kpath = np.array(kp, dtype=np.float64)

            if re.match(r'^kplane_bulk$', lw):
                bc = np.array(lines[i+1].split(), dtype=np.float64)
                bv = np.array([ l.split() for l in lines[i+2:i+4] ], \
                              dtype=np.float64)
                if bc.shape != (3,) or bv.shape != (2, 3):
                    print("KPLANE_BULK is wrong")
                    sys.exit(0)
                kp0 = bc-(bv[0]+bv[1])/2
                kp1 = kp0 + bv[0]
                kp2 = kp1 + bv[1]
                kp3 = kp0 + bv[1]
                kpath = np.array([kp0, kp1, kp2, kp3, kp0])
                kpath_name = ["O", "vec1", "vec1+2", "vec2", "O"]

            if re.match(r'^kcube_bulk$', lw):
                bc = np.array(lines[i+1].split(), dtype=np.float64)
                bv = np.array([ l.split() for l in lines[i+2:i+5] ], \
                              dtype=np.float64)
                if bc.shape != (3,) or bv.shape != (3, 3):
                    print("KCUBE_BULK is wrong")
                    sys.exit(0)
                kpath = cube_path(bc, bv)
                kpath_name = [ False ] * len(kpath)

        #----- cellの単位変更 -----#
        if   'angstrom'  in BZ_base: pass
        elif 'bohr'      in BZ_base: cell = bohrang(cell)
        elif 'ibrav'     in BZ_base:
            # 角度指定なし
            if ibrav in [1, 2, 3, -3, 4, 6, 7, 8, 9, -9, 91, 10, 11]:
                cell = ibravcell(ibrav, A = A, B = B, C = C, \
                                 cosAB = cosAB, cosAC = cosAC, cosBC = cosBC)
            # 角度指定あり
            # cosAB, cosAC, cosBCを使った角度指定
            elif ibrav in [5, -5, 12, -12, 14] and \
                (celldm4 == celldm5 == celldm6 == -2):
                cell = ibravcell(ibrav, A = A, B = B, C = C, \
                                 cosAB = cosAB, cosAC = cosAC, cosBC = cosBC)
            # celldm()を使った角度指定
            elif ibrav in [5, -5, 12, -12, 14] and \
                (cosAB == cosAC == cosBC == -2):
                if ibrav in [5, -5, 12] and (celldm5 == celldm6 == -2):
                    cell = ibravcell(ibrav, A = A, B = B, C = C, cosAB = celldm4)
                elif ibrav in [-12] and (celldm4 == celldm6 == -2):
                    cell = ibravcell(ibrav, A = A, B = B, C = C, cosAC = celldm5)
                elif ibrav in [14]:
                    cell = ibravcell(ibrav, A = A, B = B, C = C, \
                                cosAB = celldm6, cosAC = celldm5, cosBC = celldm4)
                else: _ibravError()
            elif ibrav in [13, -13]:
                #ibrav=13は計算が面倒なのでまだ非対応
                cell = ibravcell(ibrav)
            else:
                print("bz_in: ibrav not found or not support number")
                sys.exit(0)

        #----- kcellの計算 -----#
        kcell = 2*pi*(np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        if len(kpath) != 0 : kpath = np.matmul(kpath, kcell)

        #----- atom_fracの単位変更 -----#
        if   "direct"    in atom_base: pass
        elif "cartesian" in atom_base: 
            atom_frac = np.matmul( atom_frac, np.linalg.inv(cell) )

        return cell, kcell, kpath, kpath_name, atom_frac, atom_name



