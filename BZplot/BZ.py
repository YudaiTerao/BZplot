
"""
Usage:
  BZ.py [<file>] [-p|--kpath] [-v|--lcvec] [-a|--aspect]

Options:
    <file>          file名
    -p,--nokpath      kpathをplotしない
    -v,--nolcvec      逆格子ベクトルをplotしない
    -a,--aspect       BZをplotするときのaspectを揃える
"""

import sys
import re
import os 
from docopt import docopt
import numpy as np
from scipy.constants import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.quiver as qv
from mpl_toolkits.mplot3d import Axes3D
Colorlist=('#3288bd', '#d53e4f', '#fdae61', '#abdda4', '#f46d43', '#66c2a5', '#fee08b', '#5e4fa2', '#9e0142', '#e6f598')
def cminch(cm: float) -> float:
        return cm * 0.3937
def bohrang(bohr: float) -> float:
        #bohr: 5.29*10^-11m
        return bohr * physical_constants["Bohr radius"][0]*(10**10)

def manual_config():
    #file名指定なしで実行するとここが読み込まれる
    #cellはangstromで指定
  #  例: FCC, 格子定数1.28A
  #  lc = 1.28
  #  cell = np.array([[0.0, 0.5, 0.5],
  #                   [0.5, 0.0, 0.5],
  #                   [0.5, 0.5, 0.0]])
  #  cell = lc * .cell
    #ibrav指定も可能
    cell = ibravcell(3)
    kcell = 2*pi*(np.linalg.inv(cell).T)
    kpath_name, atom_name = [], []
    kpath, atom_frac = np.empty(0), np.empty(0)
    return cell, kcell, kpath, kpath_name, atom_frac, atom_name

def ibravcell(ibrav: int, lc=1, lc2=1, lc3=1):
    cell = []
    if ibrav == 1:
        cell = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])*lc
    elif ibrav == 2: 
        cell = np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]])*lc/2
    elif ibrav == 3: 
        cell = np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]])*lc/2
    elif ibrav == -3:
        cell = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])*lc/2
    elif ibrav == 4:
        cell = np.array([[1, 0, 0], [-1/2, np.sqrt(3)/2, 0], [0, 0, lc3/lc]])*lc
    return cell


#####################################################
#----- inputfileからcell,kcell,kpathの読み込み -----#
#####################################################

class BZ_input:
    def __init__(self, filename: str = ""):
        #self.cell: 実空間の格子ベクトル
        #self.kcell: k空間の逆格子ベクトル
        if filename == None: 
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = manual_config()
        elif re.match(r'^.*scf.in$', filename):
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self.read_nscf_in(filename)
        elif re.match(r'^.*win$', filename):
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self.read_win(filename)
        elif re.match(r'^.*wt.in$', filename) or \
             re.match(r'^.*weyl.in', filename):
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self.read_wt_in(filename)
        else : 
            print("BZplot.BZ_input: filename is wrong")
            sys.exit(0)

    def read_nscf_in(self, file_nscf_in):
        with open(file_nscf_in, 'r') as fn:
            #改行や行頭,行末スペース、行末のカンマやコロン、タブ、コメント行を消去
            lines = [ re.sub('^!.*$', '', re.sub('^\s+', '', re.sub('[\,\.]?\s+$', '', l))) for l in fn.readlines() ]
        ibrav = 100
        lc, lc2, lc3 = -1, -1, -1
        kpath_name, atom_name = [], []
        kpath, atom_frac = np.empty(0), np.empty(0)
        BZ_base, kp_method = "", ""

        for i, line in enumerate(lines):
            # lwはlineを小文字に統一し、コメント,空白を削除したもの
            lw = re.sub('\s', '', line).lower().split('!')[0]
            # ibravは必ず整数, 整数でなければ読み込まない
            if re.match(r'^ibrav=-?[0-9]+$', lw):
                ibrav=int(line.split("=")[1])
                for line2 in lines:
                    lw2 = re.sub('\s', '', line2).lower().split('!')[0]
                    #a, b, c, celldmは必ず実数, 実数でなければ読み込まない
                    if   re.match(r'^a=-?[0-9]+\.?[0-9]*$', lw2):
                        lc  = float(line2.split("=")[1])
                    elif re.match(r'^b=-?[0-9]+\.?[0-9]*$', lw2):
                        lc2 = float(line2.split("=")[1])
                    elif re.match(r'^c=-?[0-9]+\.?[0-9]*$', lw2):
                        lc3 = float(line2.split("=")[1])
                    elif re.match(r'^celldm\(1\)=-?[0-9]+\.?[0-9]*$', lw2):
                        lc  = bohrang(float(line2.split("=")[1]))
                    elif re.match(r'^celldm\(2\)=-?[0-9]+\.?[0-9]*$', lw2):
                        lc2 = bohrang(float(line2.split("=")[1]))
                    elif re.match(r'^celldm\(3\)=-?[0-9]+\.?[0-9]*$', lw2):
                        lc3 = bohrang(float(line2.split("=")[1]))

            elif re.match(r'^cell_parameters.+$', lw):
                BZ_base = re.sub(r'^cell_parameters', '', lw)
                cell = np.array([ l.split()[:3] for l in lines[i+1:i+4] ], dtype='float64')
            elif re.match(r'^nat=[0-9]+$', lw):
                atom_num = int(line.split("=")[1])
            elif re.match(r'^atomic_positions.+$', lw):
                atom_base = re.sub(r'^atomic_positions', '', lw)
                atom_name = [ l.split()[0] for l in lines[i+1:i+1+atom_num]]
                atom_frac = np.array([ l.split()[1:4] for l in lines[i+1:i+1+atom_num]], dtype='float64')
            elif re.match(r'^k_points.+$', lw):
                kp_method = re.sub(r'^k_points', '', lw)
                if 'tpiba_b' in kp_method or 'crystal_b' in kp_method:
                    kp_num = int(lines[i+1])
                    kpath = np.array([ l.split()[:3] for l in lines[i+2:i+2+kp_num] ], dtype='float64')
                    for j in range(i+2, i+2+kp_num):
                        if len(lines[j].split('!')) == 2:
                            kpath_name.append(lines[j].split('!')[1])

        #----- cellの単位変更 -----#
        #cell : 実空間の基底ベクトル,単位はangstromに統一
        #cellに関する部分は必ず出力に必要なため、errorの条件分岐を書いている
        lc_flag = 0
        if ibrav == 0:
            if 'alat' in BZ_base: 
                if lc == -1 : lc_flag = 1
                else: cell = cell * lc
            elif 'bohr'     in BZ_base: cell = bohrang(cell)
            elif 'angstrom' in BZ_base: pass
            else: 
                print("read_scf_in: CELL_PARAMETERS is wrong")
                sys.exit(0)
        elif 1 <= ibrav <= 3 or ibrav == -3: 
            if lc == -1: lc_flag = 1
            else: cell = ibravcell(ibrav, lc=lc)
        elif ibrav == 4:
            if lc == -1 or lc3 == -1: lc_flag = 1
            else: cell = ibravcell(ibrav, lc=lc, lc3=lc3)
        else:
            print("read_scf_in: ibrav not found or not support number")
            sys.exit(0)
        if lc_flag == 1:
            print("read_scf_in : lattice constant not found")
            sys.exit(0)

        #----- kcellの計算 -----#
        #kcell : K空間の基底ベクトル
        #kcellは実空間の基底ベクトルの逆行列の転置*2pi
        kcell = 2 * pi * (np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        if 'tpiba_b' in kp_method:
            kpath = kpath * 2 * pi / lc
        elif 'crystal_b' in kp_method:
            kpath = np.matmul(kpath, kcell)

        #----- atom_fracの単位変換 -----#
        if   'alat' in atom_base:
            atom_frac = np.matmul( atom_frac * lc, np.linalg.inv(cell) )
        elif 'bohr' in atom_base:
            atom_frac = np.matmul( bohrang(atom_frac), np.linalg.inv(cell) )
        elif 'angstrom' in atom_base:
            atom_frac = np.matmul( atom_frac, np.linalg.inv(cell) )
        elif 'crystal' in atom_base: pass

        return cell, kcell, kpath, kpath_name, atom_frac, atom_name


    def read_win(self, file_win):
        kpath_name, atom_name = [], []
        kpath, atom_frac = np.empty(0), np.empty(0)
        BZ_base, atom_base = "", ""
        cell_flag = 0

        with open(file_win, 'r') as fw:
            #改行や行頭,行末スペース、タブ、コメントを消去
            lines = [ re.sub('^!.*$', '', re.sub('^\s+', '', re.sub('\s+$', '', l))).split('!')[0] for l in fw.readlines() ]

        for i, line in enumerate(lines):
            # lwはlineを小文字に統一したもの
            lw = line.lower()

            if re.match(r'^begin\s+unit_cell_cart$', lw):
                BZ_base = lines[i+1]
                cell = np.array([ l.split() for l in lines[i+2:i+5] ], \
                                dtype='float64')
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
                kpath = np.array(kp, dtype='float64')
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
                atom_frac = np.array(af, dtype='float64')

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


    def read_wt_in(self, file_wt_in):
        with open(file_wt_in, 'r') as fwt:
            lines = fwt.readlines()
        kpath, kpath_name, atom_frac, atom_name = [], [], [], []
        for i, line in enumerate(lines):
            if "LATTICE" in line:
                cell = np.array([ l.split()[:3] for l in lines[i+2:i+5] ], dtype='f8')
                BZ_base = lines[i+1].split()[0].lower()
            if "ATOM_POSITIONS" in line:
                atom_num = int(lines[i+1].split()[0])
                atom_base = lines[i+2].split()[0].lower()
                atom_name = [ l.split()[0] for l in lines[i+3:i+3+atom_num] ]
                atom_frac = np.array([ l.split()[1:4] for l in lines[i+3:i+3+atom_num] ], dtype='f8')

            if "KPLANE_BULK" in line:
                bc = np.array(lines[i+1].split()[:3], dtype='f8')
                bv = np.array([ l.split()[:3] for l in lines[i+2:i+4] ], dtype='f8')
                kp0 = bc-(bv[0]+bv[1])/2
                kp1 = kp0 + bv[0]
                kp2 = kp1 + bv[1]
                kp3 = kp0 + bv[1]
                kpath = np.array([kp0, kp1, kp2, kp3, kp0])
                kpath_name = ["O", "vec1", "vec1+2", "vec2", "O"]

            if "KCUBE_BULK" in line:
                bc = np.array(lines[i+1].split()[:3], dtype='f8')
                bv = np.array([ l.split()[:3] for l in lines[i+2:i+5] ], dtype='f8')
                kpath = cube_path(bc, bv)
                kpath_name = [ False ] * len(kpath)

        #----- cellの単位変更 -----#
        if   'angstrom'  in BZ_base: pass
        elif 'bohr'      in BZ_base: cell = bohrang(cell)

        #----- kcellの計算 -----#
        kcell = 2*pi*(np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        if len(kpath) != 0 : kpath = np.matmul(kpath, kcell)

        #----- atom_fracの単位変更 -----#
        if   "direct"    in atom_base: pass
        elif "cartesian" in atom_base: 
            atom_frac = np.matmul( atom_frac, np.linalg.inv(cell) )

        return cell, kcell, kpath, kpath_name, atom_frac, atom_name


    def kakunin(self):
    #kcellとcellの確認用
    #素直にcellから外積や内積で計算した逆格子とkcellを比較する
        print(self.cell)
        b1=np.cross(self.cell[1], self.cell[2])
        b2=np.cross(self.cell[2], self.cell[0])
        b3=np.cross(self.cell[0], self.cell[1])
        AA=np.dot(self.cell[0], np.cross(self.cell[1], self.cell[2]))

        print((b1/AA)*2*pi)
        print((b2/AA)*2*pi)
        print((b3/AA)*2*pi)
        print(self.kcell)


#############################################
#----- 与えられたcellのvoronoi図の作成 -----#
#############################################

def get_brillouin_zone_3d(cell):
    """
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points. 

    https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell
    https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#voronoi-diagrams
    """
    cell = np.asarray(cell, dtype=float)
    assert cell.shape == (3, 3)

    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

    from scipy.spatial import Voronoi
    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []

    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        # WHY 13 ????
        # The Voronoi ridges/facets are perpendicular to the lines drawn
        # between the input points. The 14th input point is [0, 0, 0].

        if(pid[0] == 13 or pid[1] == 13):
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(vor.vertices[rid])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))

    return vor.vertices[bz_vertices], bz_ridges, bz_facets



########################################
#----- 3次元でのkpath, BZ等のplot -----#
########################################

def BZ_plot(ax, cell: np.ndarray, aspect=False):
#get_brillouinzone_3dは3つのベクトルに対してvoronoi図を作成する
#このためcellとして格子ベクトルを与えると実空間のBZを、逆格子ベクトルを与えるとk空間のBZをplotする。
    v, e, f = get_brillouin_zone_3d(cell)
    for bzp in e:
    #xx[:,0]はe[[a1,a2,a3], [b1,b2,b3], ...]の一列目(a1,b1,c1...)をすべて取り出すという意味
        ax.plot(bzp[:, 0], bzp[:, 1], bzp[:, 2], color='black', lw=1.0)

    if aspect == False:
        b1, b2, b3 = np.linalg.norm(cell, axis=1)   
        ax.set_xlim(-b1, b1)
        ax.set_ylim(-b2, b2)
        ax.set_zlim(-b3, b3)
        ax.set_box_aspect((1,1,1))
    else:
        amax = []
        for i, axis in enumerate(['x', 'y', 'z']):
            amax.append(max(abs(cell[:, i])))
            eval("ax.set_{}lim".format(axis))(-1*amax[i], amax[i])
        ax.set_box_aspect((1, amax[1]/amax[0], amax[2]/amax[0]))

def cell_plot(ax, fig, cell: np.ndarray, atom_frac, atom_name):
    clpath = cube_path(0.0, cell)
    ax.plot(clpath[:, 0], clpath[:, 1], clpath[:, 2], color='black', lw=1.0)
    amax, amin = [], []
    for i, axis in enumerate(['x', 'y', 'z']):
        amx = max(max(clpath[:, i]), 0)
        amn = min(min(clpath[:, i]), 0)
        amax.append(amx+(amx-amn)*0.4)
        amin.append(amn-(amx-amn)*0.4)
        eval("ax.set_{}lim".format(axis))(amin[i], amax[i])
    ax.set_box_aspect((1, (amax[1]-amin[1])/(amax[0]-amin[0]), (amax[2]-amin[2])/(amax[0]-amin[0])))

    ###### 各原子の色分け #####
    i, atom_dict = 0, {}
    for an in atom_name:
        if an not in atom_dict:
            atom_dict[an] = Colorlist[i]
            fig.text(0.05, 0.1+i*0.05, an, color = Colorlist[i], fontsize=20)
            i += 1
    ###############

    ##### 与えられた原子座標から  ##########################
    ##### 実空間で等価な、最も原点に近い原子座標を取得 #####
    for i, atf in enumerate(atom_frac):
        for j, af in enumerate(atf):
            if 0.9997 > af >= -0.0005: continue
            else:
                new_af = af
                if   af >= 0.9997:
                    while new_af >= 0.9997: new_af -= 1.0000
                elif af < -0.0005:
                    while new_af < -0.0005: new_af += 1.0000
                atom_frac[i][j]=new_af
    ###############

    ##### 原子の複製 #####
    natom_frac = atom_frac.copy()
    natom_name = atom_name.copy()

    for i, atf in enumerate(atom_frac):
        #new_atfには複製された原子座標が追加される。
        new_atf = [ atf.copy() ]

        for j, af in enumerate(atf):
        #複製は各軸ごとに行う。
        #例えばx軸方向の複製では,new_atf内の原子のx座標が0.0の時
        #格子ベクトル1つ分足したx座標が1.0の位置にも原子はあるはずである
        #このx座標が1.0の原子をnew_atfに追加し、
        #新しいnew_atfでy軸についても同じことを行う
            add_atf = []
            for naf in new_atf :
                if 0.0005 >= af >= -0.0005:
                    cp_naf = naf.copy()
                    cp_naf[j] += 1.0
                    add_atf.append(cp_naf)
                    natom_frac = np.append(natom_frac, cp_naf)
                    natom_name.append(atom_name[i])
            new_atf = new_atf + add_atf

    atom_coord = np.matmul(natom_frac.reshape([-1,3]), cell)
    ###############

    ##### 原子のplot #####
    radius = max([ amax[i]-amin[i] for i in range(3) ])/30
    for ac, an in zip(atom_coord, natom_name):
        ballplot(ax, radius, ac[0], ac[1], ac[2], atom_dict[an])
        print("{0}:\t[ {1[0]:.6f}, {1[1]:.6f}, {1[2]:.6f} ]".format(an, ac))
    ###############

def lcvec_plot(ax, cell, lbl_center=False):
    for i,vec in enumerate(cell):
        quiv = ax.quiver(0, 0, 0, vec[0], vec[1], vec[2], \
                        arrow_length_ratio=0.1, lw=1.2)
        quiv.set_color(Colorlist[i*2])
        if lbl_center == True:
            ax.text(vec[0]/2, vec[1]/2, vec[2]/2, \
            "   b{0},[{1[0]:.2f},{1[1]:.2f},{1[2]:.2f}]".format(i+1,vec))
        elif lbl_center == False:
            ax.text(vec[0], vec[1], vec[2], \
            "   b{0},[{1[0]:.2f},{1[1]:.2f},{1[2]:.2f}]".format(i+1,vec))


def kpath_plot(ax, kpath: np.ndarray, kpath_name=[]):
    ax.plot(kpath[:, 0], kpath[:, 1], kpath[:, 2], color='red', lw=1.2)
    for i, kp in enumerate(kpath):
        if len(kpath_name) > i:
            if kpath_name[i] == False: continue
            ax.text(kp[0], kp[1], kp[2], kpath_name[i].replace("G","$\Gamma$"), \
                    c='black', va='bottom', zorder=100)
        else: ax.text(kp[0], kp[1], kp[2], "k{}".format(i), \
                      c='black', va='bottom', zorder=100)

def cube_path(origin: np.ndarray, base_vec: np.ndarray):
    #三次元の六面体をplotするためのkpathを出力する。
    #計算エリアのplot等に使う。

    #原点(origin)と最も離れた対角にある点(X)を往復するような形で
    #六面体の一筆書きを考える。
    #通る道に重なりが多くあるため無駄は多い。
    #(8辺のplotをするために6往復×3本=18本plotすることになる)
    OtoX = [[0, 1, 2], [1, 0, 2], [2, 0, 1]]
    XtoO = [[1, 2, 0], [0, 2, 1], [0, 1, 2]]
    kpath = np.zeros((19, 3))
    kpath[0] += origin
    i = 1
    for otx, xto in zip(OtoX, XtoO):
        for ox in otx:
            kpath[i] = kpath[i-1] + base_vec[ox]
            i += 1
        for xo in xto:
            kpath[i] = kpath[i-1] - base_vec[xo]
            i += 1
    return kpath

def ballplot(ax, r, ox, oy, oz, color):
    u = np.linspace(0, 2 * np.pi, 15)
    v = np.linspace(0, np.pi, 15)
    x = r * np.outer(np.cos(u), np.sin(v))+ox
    y = r * np.outer(np.sin(u), np.sin(v))+oy
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))+oz

    # Plot the surface
    ax.plot_surface(x, y, z,color=color,rcount=50, ccount=50, alpha=1, antialiased=False)

def is_under_ssh_connection():
    # The environment variable `SSH_CONNECTION` exists only in the SSH session.
    # https://qiita.com/take_me/items/f91a3ffcee4c103a9a03
    return 'SSH_CONNECTION' in os.environ.keys()

class plotoption():
    def __init__(self):
        args = docopt(__doc__)
        self.fn = args['<file>']
        self.kp_bool = args['--nokpath']
        self.lv_bool = args['--nolcvec']
        self.bzaspect = args['--aspect']

        if is_under_ssh_connection(): 
            mpl.use('TkAgg')
        #----- Default Font in local -----#
        plt.rcParams["font.serif"] = ['Times New Roman'] + plt.rcParams["font.serif"]
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["mathtext.fontset"] = "cm"   #texfont
        plt.rcParams['font.size'] = 12

def BZplot_main():
    op = plotoption()
    bz = BZ_input(filename=op.fn)
    #--- figとaxesの作成 ---#
    fig = plt.figure(figsize=(cminch(20),cminch(20)))
    ax = fig.add_axes([ 0.05, 0.05, 0.9, 0.9], projection='3d')

    #--- BZのplot ---#
    BZ_plot(ax, bz.kcell, aspect=op.bzaspect)
    print("lattice vectors:")
    for cl in bz.cell.tolist():
        print("[ {0[0]:.6f}, {0[1]:.6f}, {0[2]:.6f} ]".format(cl))
    print("\nreciprocal lattice vectors:")
    for kl in bz.kcell.tolist():
        print("[ {0[0]:.6f}, {0[1]:.6f}, {0[2]:.6f} ]".format(kl))

    #--- 逆格子ベクトルのplot ---#
    if op.lv_bool != True:
        lcvec_plot(ax, bz.kcell)

    #--- kpathのplot ---#
    if len(bz.kpath) != 0 and op.kp_bool != True: 
        kpath_plot(ax, bz.kpath, bz.kpath_name)
        print("\nkpath:")
        for i, kp in enumerate(bz.kpath.tolist()):
            if len(bz.kpath_name) > i: 
                if bz.kpath_name[i] == False : continue
                kn = bz.kpath_name[i]
            else:  kn = "k{}".format(i)
            print("{0}:\t[ {1[0]:.6f}, {1[1]:.6f}, {1[2]:.6f} ]".format(kn, kp))
    plt.show()

def cellplot_main():
    op = plotoption()
    bz = BZ_input(filename=op.fn)
    #--- figとaxesの作成 ---#
    fig = plt.figure(figsize=(cminch(20),cminch(20)))
    ax = fig.add_axes([ 0.05, 0.05, 0.9, 0.9], projection='3d')

    cell_plot(ax, fig, bz.cell, bz.atom_frac, bz.atom_name)
    print("\nlattice vectors:")
    for cl in bz.cell.tolist():
        print("[ {0[0]:.6f}, {0[1]:.6f}, {0[2]:.6f} ]".format(cl))

    #--- 格子ベクトルのplot ---#
    if op.lv_bool != True:
        lcvec_plot(ax, bz.cell, lbl_center=True)

    plt.show()

if __name__ == '__main__' : BZplot_main()

