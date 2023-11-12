
"""
Usage:
  BZ.py [<file>] [-p|--kpath] [-v|--lcvec]

Options:
    <file>          file名
    -p,--nokpath      kpathをplotしない
    -v,--nolcvec      逆格子ベクトルをplotしない
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
Color_list=('#3288bd', '#d53e4f', '#fdae61', '#abdda4', '#f46d43', '#66c2a5', '#fee08b', '#5e4fa2', '#9e0142', '#e6f598')
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
    kpath=np.empty(0)
    kpath_name=[]
    return cell, kcell, kpath, kpath_name

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
    def __init__(self, filename=""):
        #self.cell: 実空間の格子ベクトル
        #self.kcell: k空間の逆格子ベクトル
        if filename == None: 
            self.cell, self.kcell, self.kpath, self.kpath_name \
                                             = manual_config()
        elif ".scf.in" in filename or "nscf.in" in filename: 
            self.cell, self.kcell, self.kpath, self.kpath_name \
                                 = self.read_nscf_in(filename)
        elif ".win" in filename: 
            self.cell, self.kcell, self.kpath, self.kpath_name \
                                 = self.read_win(filename)
        else : 
            print("filename is wrong")
            sys.exit(1)

    def read_nscf_in(self, file_nscf_in):
        with open(file_nscf_in, 'r') as fn:
            lines = [ line.rstrip(" ,\n") for line in fn.readlines() ]
        cell, kcell, kpath, kpath_name = [], [], [], []
        BZ_base, kp_method = "", ""

        for i,line in enumerate(lines):
            lline = line.lower().strip().replace(" ", "") #行全体を小文字化し空白等消去
            if 'ibrav' in lline:
                ibrav=int(line.split("=")[1])
            #matchを使うことで行頭のa=のみを判別
            elif re.match('a=', lline ):
                lc = float(line.split("=")[1])
            elif 'celldm(1)' in lline:
                lc = bohrang(float(line.split("=")[1]))
            elif re.match('b=', lline ):
                lc2 = float(line.split("=")[1])
            elif 'celldm(2)' in lline:
                lc2 = bohrang(float(line.split("=")[1]))
            elif re.match('c=', lline ):
                lc3 = float(line.split("=")[1])
            elif 'celldm(3)' in lline:
                lc3 = bohrang(float(line.split("=")[1]))
            elif 'CELL_PARAMETERS' in line:
                BZ_base = lline
                for j in range(i+1,i+4):
                    cell.append([float(lines[j].split()[x]) for x in range(3)])
            elif 'K_POINTS' in line:
                kp_method = lline
                if 'tpiba_b' in lline or 'crystal_b' in lline:
                    for j in range(i+2, i+2+int(lines[i+1])):
                        linelist=lines[j].split()
                        kpath.append([ float(linelist[x]) for x in range(3) ])
                        if len(linelist) > 4:
                            kpath_name.append(lines[j].replace("!","").split()[-1])

        #----- cellの単位変更 -----#
        #cell : 実空間の基底ベクトル,単位はangstromに統一
        if ibrav == 0:
            if   'alat' in BZ_base: cell = np.array(cell)*lc
            elif 'bohr' in BZ_base: cell = bohrang(np.array(cell))
            elif 'angstrom' in BZ_base: cell = np.array(cell)
        elif 1 <= ibrav <= 3 or ibrav == -3: cell = ibravcell(ibrav, lc=lc)
        elif ibrav == 4: cell = ibravcell(ibrav, lc=lc, lc3=lc3)

        #----- kcellの計算 -----#
        #kcell : K空間の基底ベクトル
        #kcellは実空間の基底ベクトルの逆行列の転置*2pi
        kcell = 2*pi*(np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        if 'tpiba_b' in kp_method:
            kpath = np.array(kpath)*2*pi/lc
        elif 'crystal_b' in kp_method:
            kpath = np.matmul(np.array(kpath), kcell)
        return cell, kcell, kpath, kpath_name

    def read_win(self, file_win):
        with open(file_win, 'r') as fw:
            lines = fw.readlines()
        cell, kcell, kpath, kpath_name = [], [], [], []

        for i, line in enumerate(lines):
            lline = line.lower()
            if 'begin unit_cell_cart' in lline:
                BZ_base = lines[i+1]
                for j in range(i+2, i+5):
                    cell.append([float(lines[j].split()[x]) for x in range(3)])
            if 'begin kpoint_path' in lline:
                j = 1
                while 'end kpoint_path' not in lines[i+j].lower():
                    if j == 1:
                        kpath_name.append(lines[i+j].split()[0])
                        kpath.append([float(lines[i+j].split()[x]) for x in range(1,4)])
                    kpath_name.append(lines[i+j].split()[-4])
                    kpath.append([float(lines[i+j].split()[x]) for x in range(-3,0)])
                    j = j + 1
        cell = np.array(cell)

        #----- kcellの計算 -----#
        kcell = 2*pi*(np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        kpath = np.matmul(np.array(kpath), kcell)

        return cell, kcell, kpath, kpath_name

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

def BZ_plot(ax, cell: np.ndarray):
#get_brillouinzone_3dは3つのベクトルに対してvoronoi図を作成する
#このためcellとして格子ベクトルを与えると実空間のBZを、逆格子ベクトルを与えるとk空間のBZをplotする。
    v, e, f = get_brillouin_zone_3d(cell)
    for bzp in e:
    #xx[:,0]はe[[a1,a2,a3], [b1,b2,b3], ...]の一列目(a1,b1,c1...)をすべて取り出すという意味
        ax.plot(bzp[:, 0], bzp[:, 1], bzp[:, 2], color='black', lw=1.0)

    b1, b2, b3 = np.linalg.norm(cell, axis=1)   
    ax.set_xlim(-b1, b1)
    ax.set_ylim(-b2, b2)
    ax.set_zlim(-b3, b3)

def lcvec_plot(ax, cell):
    for i,vec in enumerate(cell):
        quiv = ax.quiver(0, 0, 0, vec[0], vec[1], vec[2], \
                        arrow_length_ratio=0.05, lw=1)
        quiv.set_color(Color_list[i*2])
        ax.text(vec[0], vec[1], vec[2], \
        "b{0},[{1[0]:.2f},{1[1]:.2f},{1[2]:.2f}]".format(i+1,vec))

def kpath_plot(ax, kpath: np.ndarray, kpath_name=[]):
    for i in range(len(kpath)):
        if i==0: continue
        ax.quiver(kpath[i-1][0], kpath[i-1][1], kpath[i-1][2], \
                  kpath[i][0]-kpath[i-1][0], kpath[i][1]-kpath[i-1][1], \
                  kpath[i][2]-kpath[i-1][2], arrow_length_ratio=0, color="red", lw=1.8)
    for i,kp in enumerate(kpath):
        if len(kpath_name) > i:
            ax.text(kp[0], kp[1], kp[2], kpath_name[i].replace("G","$\Gamma$"), \
                    c='black', va='bottom', zorder=100)
        else: ax.text(kp[0], kp[1], kp[2], "k{}".format(i), \
                      c='black', va='bottom', zorder=100)

def is_under_ssh_connection():
    # The environment variable `SSH_CONNECTION` exists only in the SSH session.
    # https://qiita.com/take_me/items/f91a3ffcee4c103a9a03
    return 'SSH_CONNECTION' in os.environ.keys()

def main():
    args = docopt(__doc__)

    if is_under_ssh_connection(): 
        mpl.use('TkAgg')
        #----- Default Font in remote -----#
        plt.rcParams["mathtext.fontset"] = "cm"   #texfont
        plt.rcParams['font.size']=12
    else :
        #----- Default Font in local -----#
        plt.rcParams["font.serif"] = "Times New Roman"
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["mathtext.fontset"] = "cm"   #texfont
        plt.rcParams['font.size']=12

    #--- figとaxesの作成 ---#
    fig = plt.figure(figsize=(cminch(20),cminch(18)))
    ax = fig.add_axes([ 0.05, 0.05, 0.9, 0.9], projection='3d')
    ax.set_aspect('equal')

    #--- BZのplot ---#
    bz = BZ_input(filename=args['<file>'])
    BZ_plot(ax, bz.kcell)
    print("lattice vectors:")
    for cl in bz.cell.tolist():
        print("[ {0[0]:.6f}, {0[1]:.6f}, {0[2]:.6f} ]".format(cl))
    print("\nreciprocal lattice vectors:")
    for kl in bz.kcell.tolist():
        print("[ {0[0]:.6f}, {0[1]:.6f}, {0[2]:.6f} ]".format(kl))

    #--- 逆格子ベクトルのplot ---#
    if args['--nolcvec'] != True:
        lcvec_plot(ax, bz.kcell)

    #--- kpathのplot ---#
    if len(bz.kpath) != 0 and args['--nokpath'] != True: 
        kpath_plot(ax, bz.kpath, bz.kpath_name)
        print("\nkpath:")
        for i, kp in enumerate(bz.kpath.tolist()):
            if len(bz.kpath_name) > i: kn = bz.kpath_name[i]
            else:  kn = "k{}".format(i)
            print("{0}:\t[ {1[0]:.6f}, {1[1]:.6f}, {1[2]:.6f} ]".format(kn, kp))
    plt.show()

