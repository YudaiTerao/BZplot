
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
    kpath=np.empty(0)
    kpath_name=[]
    atom_name=[]
    atom_frac=[]
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
    def __init__(self, filename=""):
        #self.cell: 実空間の格子ベクトル
        #self.kcell: k空間の逆格子ベクトル
        if filename == None: 
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = manual_config()
        elif ".scf.in" in filename or ".nscf.in" in filename:
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self.read_nscf_in(filename)
        elif ".win" in filename or "_win" in filename:
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self.read_win(filename)
        elif "wt.in" in filename:
            self.cell, self.kcell, self.kpath, self.kpath_name, \
            self.atom_frac, self.atom_name = self.read_wt_in(filename)
        else : 
            print("filename is wrong")
            sys.exit(1)

    def read_nscf_in(self, file_nscf_in):
        with open(file_nscf_in, 'r') as fn:
            lines = [ line.rstrip(" ,\n") for line in fn.readlines() ]
        kpath, kpath_name, atom_frac, atom_name = [], [], [], []
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
                cell = np.array([ l.split()[:3] for l in lines[i+1:i+4] ], dtype='f8')
            elif re.match('nat=', lline ):
                atom_num = int(line.split("=")[1])
            elif 'ATOMIC_POSITIONS' in line:
                atom_base = lline
                atom_name = [ l.split()[0] for l in lines[i+1:i+1+atom_num]]
                atom_frac = np.array([ l.split()[1:4] for l in lines[i+1:i+1+atom_num]], dtype='f8')
            elif 'K_POINTS' in line:
                kp_method = lline
                if 'tpiba_b' in lline or 'crystal_b' in lline:
                    kp_num = int(lines[i+1])
                    kpath = np.array([ l.split()[:3] for l in lines[i+2:i+2+kp_num] ], dtype='f8')
                    for j in range(i+2, i+2+kp_num):
                        if len(line.split()) > 4:
                            kpath_name.append(lines[j].replace("!"," ").split()[-1])

        #----- cellの単位変更 -----#
        #cell : 実空間の基底ベクトル,単位はangstromに統一
        if ibrav == 0:
            if   'alat' in BZ_base: cell = cell*lc
            elif 'bohr' in BZ_base: cell = bohrang(cell)
            elif 'angstrom' in BZ_base: pass
        elif 1 <= ibrav <= 3 or ibrav == -3: cell = ibravcell(ibrav, lc=lc)
        elif ibrav == 4: cell = ibravcell(ibrav, lc=lc, lc3=lc3)

        #----- kcellの計算 -----#
        #kcell : K空間の基底ベクトル
        #kcellは実空間の基底ベクトルの逆行列の転置*2pi
        kcell = 2*pi*(np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        if 'tpiba_b' in kp_method:
            kpath = kpath*2*pi/lc
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
        with open(file_win, 'r') as fw:
            lines = fw.readlines()
        kpath, kpath_name, atom_frac, atom_name = [], [], [], []

        for i, line in enumerate(lines):
            lline = line.lower()
            if 'begin unit_cell_cart' in lline:
                BZ_base = lines[i+1]
                cell = np.array([ l.split()[:4] for l in lines[i+2:i+5] ], dtype='f8')
            elif 'begin kpoint_path' in lline:
                j = 1
                while 'end kpoint_path' not in lines[i+j].lower():
                    kpath_name.append(lines[i+j].split()[0])
                    kpath.append(lines[i+j].split()[1:4])
                    j = j + 1
                kpath_name.append(lines[i+j-1].split()[-4])
                kpath.append(lines[i+j-1].split()[-3:])
                kpath = np.array(kpath, dtype='f8')
            elif 'begin atoms_frac' in lline or 'begin atoms_cart' in lline:
                atom_base = lline
                j = 1
                while 'end atoms_' not in lines[i+j].lower():
                    atom_name.append(lines[i+j].split()[0])
                    atom_frac.append(lines[i+j].split()[1:4])
                    j = j + 1
                atom_frac = np.array(atom_frac, dtype='f8')

        #----- cellの単位変更 -----#
        if   'ang'  in BZ_base: pass
        elif 'bohr' in BZ_base: cell = bohrang(cell)

        #----- kcellの計算 -----#
        kcell = 2*pi*(np.linalg.inv(cell).T)

        #----- kpathの単位変更 -----#
        if len(kpath) != 0 : kpath = np.matmul(kpath, kcell)

        #----- atom_fracの単位変更 -----#
        if   'atoms_frac' in atom_base : pass
        elif 'atoms_cart' in atom_base : 
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


    i, atom_dict = 0, {}
    for an in atom_name:
        if an not in atom_dict:
            atom_dict[an] = Colorlist[i]
            fig.text(0.05, 0.1+i*0.05, an, color = Colorlist[i], fontsize=20)
            i += 1

    for i, atf in enumerate(atom_frac):
        for j, af in enumerate(atf):
            if 0.9999 > af >= -0.0001: continue
            else:
                if   af > 0.9999:
                    new_af = af
                    while new_af >= -0.0001 : new_af -= 1.0
                elif af < - 0.0001: new_af = af
                while new_af < -0.0001 : new_af += 1.0
                atom_frac[i][j]=new_af

    natom_frac = atom_frac.copy()
    natom_name = atom_name.copy()

    for i, atf in enumerate(atom_frac):
        new_atf = [ atf.copy() ]
        for j, af in enumerate(atf):
            newadd = []
            for naf in new_atf :
                cp_naf = naf.copy()
                if 0.0001 >= af >= -0.0001:
                    cp_naf[j] += 1.0
                    newadd.append(cp_naf)
                    natom_frac = np.append(natom_frac, cp_naf)
                    natom_name.append(atom_name[i])
            for nad in newadd : new_atf.append(nad)

    atom_coord = np.matmul(natom_frac.reshape([-1,3]), cell)

    radius = max([ amax[i]-amin[i] for i in range(3) ])/30
    for ac, an in zip(atom_coord, natom_name):
        ballplot(ax, radius, ac[0], ac[1], ac[2], atom_dict[an])
        print("{0}:\t[ {1[0]:.6f}, {1[1]:.6f}, {1[2]:.6f} ]".format(an, ac))

def lcvec_plot(ax, cell):
    for i,vec in enumerate(cell):
        quiv = ax.quiver(0, 0, 0, vec[0], vec[1], vec[2], \
                        arrow_length_ratio=0.1, lw=1.2)
        quiv.set_color(Colorlist[i*2])
        ax.text(vec[0]/2, vec[1]/2, vec[2]/2, \
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
            #----- Default Font in remote -----#
            plt.rcParams["mathtext.fontset"] = "cm"   #texfont
            plt.rcParams['font.size']=12
        else :
            #----- Default Font in local -----#
            plt.rcParams["font.serif"] = "Times New Roman"
            plt.rcParams["font.family"] = "serif"
            plt.rcParams["mathtext.fontset"] = "cm"   #texfont
            plt.rcParams['font.size']=12

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
        lcvec_plot(ax, bz.cell)

    plt.show()

if __name__ == '__main__' : BZplot_main()

