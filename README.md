# BZplot
quantumESPRESSOやWannier90, WannierToolsのinputからBrilluan zoneや結晶構造を三次元でプロットする<br><br>
<img src="https://github.com/user-attachments/assets/2631d5be-d202-4b08-8f7b-25b901cf6475" width="47.2%" align="top" />
<img src="https://github.com/user-attachments/assets/9e2822b5-6e30-48be-95c8-6a76968e7347" width="45%"  align="top"/>
## Install
  installしたい場所で
  ```
  git clone https://github.com/YudaiTerao/BZplot.git
  pip install ./BZplot
  ```
  と実行すればbz,cellコマンドが使えるようになります。<br><br>
  uninstallは
  ```
  pip uninstall BZplot
  ```
## 仕様
  ### -- Command --
  ```
  bz Fe.scf.in
  bz Fe.nscf.in
  bz Fe.win
  bz wt.in
  ```
  のコマンドでBurilluanZoneをmatplotlibで出力します。<br><br>
  ```
  cell Fe.scf.in
  cell Fe.nscf.in
  cell Fe.win
  cell wt.in
  ```
  のコマンドで実空間のunit cellをmatplotlibで出力します。
  <br><br><br>
  ### -- Options --
  -h&ensp;  :&ensp;&ensp;  オプション情報の出力<br>
  -p|--nokpath&ensp; :&ensp;&ensp;  kpathを表示しない<br>
  -v|--nolcvec&ensp; :&ensp;&ensp;  逆格子を表示しない<br>
  -b|--nobackground&ensp;  :&ensp;&ensp;  gridなどの背景を表示しない<br>
  -t|--vectext&ensp;  :&ensp;&ensp;  ベクトルのラベルや座標を表示しない<br>
    <br><br><br>
### -- Tips --
- `nscf.in`のkpathのラベルは下のように&ensp;!&ensp;以下に書かれたものが反映されます。
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
- inputなしでBZやkpathを出力するには、コード内の関数の`manual config`に結晶構造とkpathを直打ちし、<br>
ファイル名を引数に付けずに実行してください。<br>
<br>

## ssh先での表示
defaultでssh先かローカルかを判定し、ssh先ならばバックエンドをtkaggに変更して出力するようになっています。<br>
このため、
```
ssh -Y remote名
```
でssh接続すればssh先でも使用できます。<br>

