# BZplot
quantumESPRESSOやWannier90, WannierToolsのinputからBrilluan zoneや結晶構造を三次元でプロットする<br><br>
<img src="https://github.com/user-attachments/assets/2631d5be-d202-4b08-8f7b-25b901cf6475" width="37.2%" align="top" />
<img src="https://github.com/user-attachments/assets/9e2822b5-6e30-48be-95c8-6a76968e7347" width="35.5%"  align="top"/>
# Install
  インストールしたい場所で以下のコマンドを打つことでインストールできる
  ```
  git clone https://github.com/YudaiTerao/BZplot.git
  pip install ./BZplot
  ```
  必要なパッケージはsetup.pyのコメント内に書いてあるが、自動的にダウンロードしないようにしている。<br>
  (anaconda環境などとの競合を避けるため。コメントを外せば自動的に必要なパッケージがpipでダウンロードされる(非推奨))<br>
  このため環境に合わせて足りないパッケージは導入する必要がある。<br><br>
  下記のどちらかのコマンドを実行して、Brilloan zoneまたは結晶格子が表示されればインストールは成功である。
  ```
  bz
  cell
  ```
  
  アンインストールは
  ```
  pip uninstall BZplot
  ```
  
# 仕様
  ## Functions
  具体的な引数や戻り値は、コード内のコメントを参照<br>
  - `BZplot.ibarvcell`<br>
  &ensp;格子定数とibrav番号を与えると、基本並進ベクトルを戻り値で返す。ibrav=13,-13以外に対応<br><br>
  - `BZplot.bz_plot`<br>
  &ensp;三本の一次独立な基底ベクトルを与えると、voronoi図を求め、 Brilluan zoneをプロットする<br><br>
  - `BZplot.get_brillouin_zone_3d`<br>
  &ensp;三本の基底ベクトルを与えると、voronoi図を求める<br><br>
  - `BZplot.lcvec_plot`<br>
  &ensp;格子ベクトル, 逆格子ベクトルをプロットする<br><br>
  - `BZplot.kpath_plot`<br>
  &ensp;kpathをプロットする<br><br>
  - `BZplot.cell_plot`<br>
  &ensp;三本の基底ベクトルを与えると、それに囲まれた六面体をプロットする<br><br>
  - `BZplot.cube_path`<br>
  &ensp;六面体を一筆書きするような経路を出力する<br><br>
  - `BZplot.atom_copy`<br>
  &ensp;与えられた原子位置から、それと等価な原子位置を求め、複製する<br><br>
  - `BZplot.atom_plot`<br>
  &ensp;原子をプロットする<br><br>
  - `BZplot.ballplot`<br>
  &ensp;球体をプロットする<br><br>
  ## Command
  インストールすると自動的に`bz`, `cell`コマンドが使えるようになる。<br>
  これらは引数にQuantumESPRESSOやWannier90, Wannier toolsのinput fileを与えると、拡張子からfileを判別して読みこみ、  Brilluan zoneや結晶構造, kpath, 計算領域を出力する。<br>
  実行コマンドの例は以下の通り.
  ```
  bz Fe.scf.in　　
  bz wt.in
  cell Fe.nscf.in
  cell Fe.win
  ```
  具体的には次のように判別される。<br>
  &ensp;&ensp;`.scf.in`, `.nscf.in` :&ensp;QuantumESPRESSOのinputと判別<br>
  &ensp;&ensp;`.win` :&ensp;Wannier90のinputと判別<br>
  &ensp;&ensp;`.in` :&ensp;WannierToolsのinputと判別<br>
　<br>
  ### Options
  コマンドのオプションは以下の通り。<br><br>
  &ensp;&ensp;`-h`&ensp;  :&ensp;&ensp;  オプション情報の出力<br>
  &ensp;&ensp;`-p|--nokpath`&ensp; :&ensp;&ensp;  kpathを表示しない<br>
  &ensp;&ensp;`-v|--nolcvec`&ensp; :&ensp;&ensp;  逆格子を表示しない<br>
  &ensp;&ensp;`-b|--nobackground`&ensp;  :&ensp;&ensp;  gridなどの背景を表示しない<br>
  &ensp;&ensp;`-t|--vectext`&ensp;  :&ensp;&ensp;  ベクトルのラベルや座標を表示しない<br>
  
  &ensp;&ensp;`-x　<xbd>`, `-y <ybd>`, `-z <zbd>`&ensp;  :<br>
  &ensp;&ensp;&ensp;&ensp;`cell`専用オプション, 各軸方向に表示する単位胞の数を指定する。<br>
  &ensp;&ensp;&ensp;&ensp;カンマで区切り、最小値と最大値を指定する。<br>
  &ensp;&ensp;&ensp;&ensp;defaultはすべて0,1で、この場合単位胞は一つだけ出力される。
  <br><br>
  実行コマンドの例は以下の通り.
  ```
  bz Fe.nscf.in -p -v -b
  cell Fe.nscf.in -x 0,2 -y 0,2 -z 0,1
  ```
  
  ### Tips
  - `.nscf.in`のkpathのラベルは下のように&ensp;!&ensp;以下に書かれたものが反映される。

```
K_POINTS {crystal_b}
8
    0.0000000000  0.0000000000  0.0000000000    20  !  G
    0.5000000000  0.0000000000  0.0000000000    20  !  M
    0.3333333333  0.3333333333  0.0000000000    20  !  K
    0.0000000000  0.0000000000  0.0000000000    20  !  G
    0.0000000000  0.0000000000  0.5000000000    20  !  A
    0.5000000000  0.0000000000  0.5000000000    20  !  L
    0.3333333333  0.3333333333  0.5000000000    20  !  H
    0.0000000000  0.0000000000  0.5000000000    20  !  A
```

  ### ssh先でのコマンドによる表示
  defaultでssh先かローカルかを判定し、ssh先ならばバックエンドをtkaggに変更して出力するようになっている。<br>
  このため、
```
ssh -Y remote名
```
  でssh接続すればssh先でも使用できる。<br>

