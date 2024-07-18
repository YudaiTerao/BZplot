
from typing import List
import numpy as np
from matplotlib.figure import Figure
import matplotlib.quiver as qv
from mpl_toolkits.mplot3d import Axes3D
Colorlist=('#3288bd', '#d53e4f', '#fdae61', '#abdda4', '#f46d43', '#66c2a5', '#fee08b', '#5e4fa2', '#9e0142', '#e6f598')

def adjust_aspect(ax: Axes3D):
    axis_len = np.zeros(3, dtype = np.float64)
    axis_rate = np.zeros(3, dtype = np.float64)
    for i, axis in enumerate(['x', 'y', 'z']):
        amn, amx = eval("ax.get_{}lim".format(axis))()
        axis_len[i] = amx - amn
        eval("ax.set_{}lim".format(axis))(amn-axis_len[i]*0.2, amx+axis_len[i]*0.2)
    for i, axis in enumerate(['x', 'y', 'z']):
        if i == axis_len[i].argmin():
            len_base = axis_len[i]
            axis_rate[i] = 1
        else: axis_rate[i] = axis_len[i]/len_base
    ax.set_box_aspect(tuple(axis_rate))

def bz_plot(
    ax: Axes3D,
    kcell: np.ndarray,
    bz_linewidth: float = 1.0,
    aspect: bool = True,
    background: bool = True
):
    """
    与えられた基底ベクトルから3次元のvoronoi図をプロットする

    Args:
        ax(Axes3D): プロットを記録するAxesオブジェクト
        kcell(array): 3次元の基底ベクトル(3, 3)
                      逆格子ベクトルを入れればBrilloan zoneと等価になる
        bz_linewidth(float): BZの境界線の太さ
        aspect(bool): aspect比を揃えるかどうか
        background(bool): 背景の目盛線を表示するかどうか
    """
    v, e, f = get_brillouin_zone_3d(kcell)
    for bzp in e:
        ax.plot(bzp[:, 0], bzp[:, 1], bzp[:, 2], color='black', lw=bz_linewidth)

#    if aspect == False:
#        b1, b2, b3 = np.linalg.norm(kcell, axis=1)
#        xm, ym, zm = kcell.max(axis=0)
#        ax.set_xlim(-xm, xm)
#        ax.set_ylim(-ym, ym)
#        ax.set_zlim(-zm, zm)
#        # 各軸の長さ比を揃える, 1目盛のサイズ比は等しくない
#        ax.set_box_aspect((1,1,1))
#    elif aspect == True:
#        # 各軸の長さに合わせてaspect比を調整する
#        amax = []
#        for i, axis in enumerate(['x', 'y', 'z']):
#            amax.append(max(abs(kcell[:, i])))
#            eval("ax.set_{}lim".format(axis))(-1*amax[i], amax[i])
#        ax.set_box_aspect((1, amax[1]/amax[0], amax[2]/amax[0]))
    if background == False:
        for axis in ['x', 'y', 'z']:
            eval("ax.{}axis.set_pane_color".format(axis))('white', alpha=0.0)
        ax.grid(None)

def get_brillouin_zone_3d(cell: np.ndarray):
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

def lcvec_plot(
    ax: Axes3D, 
    cell: np.ndarray,
    lbl: bool = True,
    lbl_center: bool = False
):
    """
    Axes3Dに任意個数のベクトル追加し、矢印でプロットする

    Args:
        ax(Axes3D): プロットを記録するAxesオブジェクト
        cell(array): 任意個数のベクトルを格納したArray(, 3)
        lbl_center(bool): ベクトル成分のテキストをベクトル矢印の先端に
                             表示するか(False)、ベクトル矢印の中心に表示す
                             るか(True)
    """
    for i,vec in enumerate(cell):
        quiv = ax.quiver(0, 0, 0, vec[0], vec[1], vec[2], \
                        arrow_length_ratio=0.1, lw=1.2)
        quiv.set_color(Colorlist[i*2])
        if lbl_center == True:
            if lbl == True: ax.text(vec[0]/2, vec[1]/2, vec[2]/2, \
            "   b{0},[{1[0]:.2f},{1[1]:.2f},{1[2]:.2f}]".format(i+1,vec))
        elif lbl_center == False:
            if lbl == True: ax.text(vec[0], vec[1], vec[2], \
            "   b{0},[{1[0]:.2f},{1[1]:.2f},{1[2]:.2f}]".format(i+1,vec))

def kpath_plot(
    ax: Axes3D, 
    kpath: np.ndarray,
    kpath_name: List = [],
    kp_linewidth: float = 1.2
):
    """
    kpathをプロットする

    Args:
        ax(Axes3D): プロットを記録するAxesオブジェクト
        kpath(array): kpathの座標(cart)(, 3)
        kpath_name(List): kpathのラベル(:)
        kp_linewidth: kpathの線の太さ
    """
    ax.plot(kpath[:, 0], kpath[:, 1], kpath[:, 2], color='red', lw=kp_linewidth)
    for i, kp in enumerate(kpath):
        if len(kpath_name) > i:
            if kpath_name[i] == False: continue
            ax.text(kp[0], kp[1], kp[2], kpath_name[i].replace("G","$\Gamma$"), \
                    c='black', va='bottom', zorder=200+i)
        else: ax.text(kp[0], kp[1], kp[2], "k{}".format(i), \
                      c='black', va='bottom', zorder=200+i)

def cell_plot(
    ax: Axes3D,
    cell: np.ndarray,
    boundary: np.ndarray,
    cl_linewidth: float = 1.0,
    aspect: bool = True,
    background: bool = True
):
    """
    3本のベクトルに囲まれた六面体のプロット
    結晶構造の出力に用いる

    Args:
        ax(Axes3D): プロットを記録するAxesオブジェクト
        cell(array): 3次元の実空間の基底ベクトル(3, 3)
        boundary(array): unitcellの繰り返し数,各軸ごとに最小,最大を指定(3, 2)
        cl_linewidth(float): unitcellの境界線の太さ
        aspect(bool): aspect比を揃えるかどうか
        background(bool): 背景の目盛線を表示するかどうか
    """
    ucl = cube_path(np.zeros(3), cell)
    ax.plot(ucl[:, 0], ucl[:, 1], ucl[:, 2], color='black', lw=cl_linewidth)

    bd_len = boundary[:, 1] - boundary[:, 0]
    allcell = cell * bd_len.reshape(3, 1)
    origin = (cell * boundary[:, 0].reshape(3, 1)).sum(axis=0)
    cl = cube_path(origin, allcell)

#    if aspect == False:
#        cmx = cl.max(axis=0)
#        cmn = cl.min(axis=0)
#        xmax, ymax, zmax = cl.max(axis=0)*1.1
#        xmin, ymin, zmin = cl.min(axis=0)*1.1
#        ax.set_xlim(xmin, xmax)
#        ax.set_ylim(ymin, ymax)
#        ax.set_zlim(zmin, zmax)
#        # 各軸の長さ比を揃える, 1目盛のサイズ比は等しくない
#        ax.set_box_aspect((1,1,1))
#    elif aspect == True:
#        amax, amin = [], []
#        for i, axis in enumerate(['x', 'y', 'z']):
#            amx = max(max(cl[:, i]), 0)
#            amn = min(min(cl[:, i]), 0)
#            amax.append(amx+(amx-amn)*0.4)
#            amin.append(amn-(amx-amn)*0.4)
#            eval("ax.set_{}lim".format(axis))(amin[i], amax[i])
#        ax.set_box_aspect((1, (amax[1]-amin[1])/(amax[0]-amin[0]), \
#                              (amax[2]-amin[2])/(amax[0]-amin[0])))
    if background == False:
        for axis in ['x', 'y', 'z']:
            eval("ax.{}axis.set_pane_color".format(axis))('white', alpha=0.0)
        ax.grid(None)

def cube_path(origin: np.ndarray, base_vec: np.ndarray):
    """
    三次元の六面体をプロットするためのkpathを出力する
    六面体の一筆書きをkpathとして出力する。
    BZ.pyのkpath_plotと組み合わせて、KCUBE_BULKなどの計算領域を表示できる

    原点(origin)と最も離れた対角にある点(X)を往復するような形で
    六面体の一筆書きを考えている。
    (8辺のplotをするために6往復×3本=18本plotすることになる)

    Args:
        origin(array): 六面体の原点座標(3)
        base_vec(array): 六面体の基底ベクトル(3, 3)
    Returns:
        kpath(array): 六面体の一筆書きのkpathの座標(:3)
    """
    OtoX = [[0, 1, 2], [1, 0, 2], [2, 0, 1]]
    XtoO = [[1, 2, 0], [0, 2, 1], [0, 1, 2]]
    kpath = np.zeros((19, 3), dtype=np.float64)
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

def atom_copy(
    atom_frac: np.ndarray,
    atom_name: List,
    boundary: np.ndarray
):
    """
    与えられた原子座標とBZ内で等価な原子座標を計算し、複製する

    Args:
        atom_frac(array): 分数の原子座標(, 3)
        atom_name(List): 原子名(, 3)
        boundary(array): unitcellの繰り返し数,各軸ごとに最小,最大を指定(3, 2)
    """

    # 与えられた原子座標から実空間で等価な、最も原点に近い原子座標を取得
    for i, atf in enumerate(atom_frac):
        for j, af in enumerate(atf):
            if 0.9997 > af >= -0.0005: continue
            else:
                af0 = af
                if   af >= 0.9997:
                    while af0 >= 0.9997: af0 -= 1.0000
                elif af < -0.0005:
                    while af0 < -0.0005: af0 += 1.0000
                atom_frac[i][j] = af0

    # 原子の複製
    new_atom_frac = atom_frac.copy().tolist()
    new_atom_name = atom_name.copy()

    for na, origin_atomfrac in enumerate(atom_frac):
        #new_atfには複製された原子座標が追加される。
        copy_atomsfrac = [ origin_atomfrac.tolist() ]
        temp_atomsfrac = [ origin_atomfrac.tolist() ]

        for i, af in enumerate(origin_atomfrac):
        #複製は各軸ごとに行う。
        #例えばx軸方向の複製では,new_atf内の原子のx座標が0.0の時
        #格子ベクトル1つ分足したx座標が1.0の位置にも原子はあるはずである
        #このx座標が1.0の原子をnew_atfに追加し、
        #新しいnew_atfでy軸についても同じことを行う
            for bd in range(boundary[i][0], boundary[i][1] + 1):
                if bd == 0: continue
                if bd == boundary[i][1] and \
                   not 0.0005 >= af >= -0.0005: continue
                caf = np.array(copy_atomsfrac)
                caf[:, i] = caf[:, i] + bd
                temp_atomsfrac += caf.tolist()
                new_atom_frac += caf.tolist()
                new_atom_name += [ atom_name[na] for n in range(len(caf)) ]
            copy_atomsfrac = temp_atomsfrac.copy()

    return np.array(new_atom_frac), new_atom_name

def atom_plot(
    ax: Axes3D,
    fig: Figure,
    cell: np.ndarray,
    atom_frac: np.ndarray,
    atom_name: List = []
):
    """
    原子をプロットする

    Args:
        ax(Axes3D): プロットを記録するAxesオブジェクト
        fig(Figure): Axesオブジェクトが存在するFigureオブジェクト
        cell(array): 任意個数のベクトルを格納したArray(, 3)
        atom_frac(array): 分数の原子座標(, 3)
        atom_name(List): 原子名(, 3)
    """
    i, atom_dict = 0, {}
    # 各原子を色分けし,ラベルをfigに追加
    for an in atom_name:
        if an not in atom_dict:
            atom_dict[an] = Colorlist[i]
            fig.text(0.05, 0.1+i*0.05, an, color = Colorlist[i], fontsize=20)
            i += 1
    atom_coord = np.matmul(atom_frac.reshape([-1,3]), cell)

    # 原子をプロットするため、球体の半径を決める
    maxlen = 0
    cl = cube_path(np.zeros(3), cell)
    for i in range(3):
        amx = max(max(cl[:, i]), 0)
        amn = min(min(cl[:, i]), 0)
        maxlen = max([amx - amn, maxlen])
    radius = maxlen/17

    if len(atom_name) > 75: rc, cc = 10, 10
    else: rc, cc = 50, 50
    # 原子のプロット
    for ac, an in zip(atom_coord, atom_name):
        ballplot(ax, r=radius, o=ac, color=atom_dict[an], rcount=rc, ccount=cc)
        print("{0}:\t[{1[0]:>10.6f},{1[1]:>10.6f},{1[2]:>10.6f} ]".format(an, ac))

def ballplot(
    ax: Axes3D,
    r: float = 1.0,
    o: np.ndarray = np.array([0, 0, 0]),
    color: str = 'gray',
    rcount: int = 50,
    ccount: int = 50
):
    """
    球体のプロット

    Args:
        ax(Axes3D): プロットを記録するAxesオブジェクト
        r(float): 任意個数のベクトルを格納したArray(, 3)
        o(array): 原点座標
        color(str): 球体の色
    """
    u = np.linspace(0, 2 * np.pi, 15)
    v = np.linspace(0, np.pi, 15)
    x = r * np.outer(np.cos(u), np.sin(v)) + o[0]
    y = r * np.outer(np.sin(u), np.sin(v)) + o[1]
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + o[2]

    # Plot the surface
    ax.plot_surface(x, y, z,color=color, rcount=rcount, ccount=ccount, alpha=1, antialiased=False)



