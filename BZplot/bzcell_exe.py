
"""
Usage:
 bzcell_exe.py [<file>] [-p|--nokpath] [-v|--nolcvec] [-b|--nobackground] [-t|--vectext] [-x <xbd>] [-y <ybd>] [-z <zbd>]

Options:
    <file>             file名 [default:None]
    -p,--nokpath       kpathをplotしない
    -v,--nolcvec       逆格子ベクトルをplotしない
    -b,--nobackground  gridなどの背景を表示しない
    -t,--vectext      ベクトルのラベルや座標を表示する
    -x <xbd>    x軸方向のunitcell数 [default: 0,1]
    -y <ybd>    y軸方向のunitcell数 [default: 0,1]
    -z <zbd>    z軸方向のunitcell数 [default: 0,1]
"""
import os
import re
import sys
from docopt import docopt
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BZplot import read_input as RI
from BZplot import bzplot as BP

def cminch(cm: float) -> float:
        return cm * 0.3937

def is_under_ssh_connection():
    # The environment variable `SSH_CONNECTION` exists only in the SSH session.
    # https://qiita.com/take_me/items/f91a3ffcee4c103a9a03
    return 'SSH_CONNECTION' in os.environ.keys()

def bz():
    cmd = bzplot_cmd()
    cmd.bz_cmd()

def cell():
    cmd = bzplot_cmd()
    cmd.cell_cmd()

class bzplot_cmd():
    """
    bz, cellコマンドの実行プログラム

    Attributes:
        inputpath(Path|None): inputfileのpath
        kp_bool(bool): kpathを表示するかどうか
        lv_bool(bool): 逆格子ベクトルを表示するかどうか
        bg_bool(bool): 背景を表示するかどうか
        vectext_bool(bool): ベクトルの座標を表示するかどうか
        boundary(array): 各軸方向のunitcellの数
    """

    def __init__(self):
        args = docopt(__doc__)
        if args['<file>'] is not None:
            self.inputpath = Path(args['<file>'])
        else: self.inputpath = None
        self.kp_bool = not args['--nokpath']
        self.lv_bool = not args['--nolcvec']
        self.bg_bool = not args['--nobackground']
        self.vectext_bool = args['--vectext']
        if not re.match(r'(-[1-9],[0-9]|-?0,[1-9])', args['-x']) \
           or not re.match(r'(-[1-9],[0-9]|-?0,[1-9])', args['-y']) \
           or not re.match(r'(-[1-9],[0-9]|-?0,[1-9])', args['-z']):
            print("bzplot_cmd: boundary is wrong")
            sys.exit(0)
        else:
            self.boundary = np.zeros((3, 2), dtype = np.int8)
            self.boundary[0] = np.array(args['-x'].split(','))
            self.boundary[1] = np.array(args['-y'].split(','))
            self.boundary[2] = np.array(args['-z'].split(','))

        if is_under_ssh_connection(): 
            mpl.use('TkAgg')
        #----- Default Font -----#
        plt.rcParams["font.serif"] = ['Times New Roman'] + \
                                     plt.rcParams["font.serif"]
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["mathtext.fontset"] = "cm"   #texfont
        plt.rcParams['font.size'] = 12

    def bz_cmd(self):
        bz = RI.BZ_input(filepath = self.inputpath)
        #--- figとaxesの作成 ---#
        fig = plt.figure(figsize=(cminch(20),cminch(20)))
        ax = fig.add_axes([0.05, 0.05, 0.9, 0.9], projection='3d')

        #--- BZのplot ---#
        BP.bz_plot(ax, bz.kcell, background=self.bg_bool)
        print("lattice vectors:")
        for cl in bz.cell.tolist():
            print("[{0[0]:>10.6f}, {0[1]:>10.6f}, {0[2]:>10.6f} ]".format(cl))
        print("\nreciprocal lattice vectors:")
        for kl in bz.kcell.tolist():
            print("[{0[0]:>10.6f}, {0[1]:>10.6f}, {0[2]:>10.6f} ]".format(kl))

        #--- 逆格子ベクトルのplot ---#
        if self.lv_bool == True:
            BP.lcvec_plot(ax, bz.kcell, lbl=self.vectext_bool)

        #--- kpathのplot ---#
        if len(bz.kpath) != 0 and self.kp_bool == True: 
            BP.kpath_plot(ax, bz.kpath, bz.kpath_name)
            print("\nkpath:")
            for i, kp in enumerate(bz.kpath.tolist()):
                if len(bz.kpath_name) > i: 
                    if bz.kpath_name[i] == False : continue
                    kn = bz.kpath_name[i]
                else:  kn = "k{}".format(i)
                print("{0}:\t[{1[0]:>10.6f},{1[1]:>10.6f},{1[2]:>10.6f} ]".format(kn, kp))
        BP.adjust_aspect(ax)
        plt.show()

    def cell_cmd(self):
        bz = RI.BZ_input(filepath = self.inputpath)
        #--- figとaxesの作成 ---#
        fig = plt.figure(figsize=(cminch(20),cminch(20)))
        ax = fig.add_axes([ 0.05, 0.05, 0.9, 0.9], projection='3d')
        print("lattice vectors:")
        for cl in bz.cell.tolist():
            print("[{0[0]:>10.6f},{0[1]:>10.6f}, {0[2]:>10.6f} ]".format(cl))
        print('')

        BP.cell_plot(ax, bz.cell, self.boundary, background=self.bg_bool)
        new_af, new_an = BP.atom_copy(bz.atom_frac, bz.atom_name, self.boundary)
        BP.atom_plot(ax, fig, bz.cell, new_af, new_an)

        #--- 格子ベクトルのplot ---#
        if self.lv_bool == True:
            BP.lcvec_plot(ax, bz.cell, lbl=self.vectext_bool, lbl_center=True)
        BP.adjust_aspect(ax)
        plt.show()



