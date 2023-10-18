# BZplot
quantumESPRESSOやWannier90のinputからBrilluan zoneを3Dplotし、逆格子, kpathを表示する  
![BZplot_sample](https://github.com/YudaiTerao/BZplot/assets/103988651/b9f19c1a-81b0-48b2-9ed6-98b62b49bdac)
## 仕様
  python3 BZ.py&ensp; -v&ensp; `Fe.scf.in`<br>
  python3 BZ.py&ensp; -p&ensp; -v&ensp; `Fe.nscf.in`（bandplot用のnscf.in)<br>
  python3 BZ.py&ensp; -p&ensp; -v&ensp; `Fe.win`<br>
  <br>
  -p|--kpath&ensp; :&ensp;&ensp;  kpathを表示する<br>
  -v|--vec&ensp;   :&ensp;&ensp;  逆格子を表示する<br>
  <br>
- `nscf.in`のkpathをplotするとき、kpathのラベルは下のように&ensp;!&ensp;以下に書かれたものが反映されます。<br><br>
K_POINTS {crystal_b}<br>
8<br>
    &emsp;0.0000000000&emsp;0.0000000000&emsp;0.0000000000&ensp;  20&ensp;  !&ensp;  G<br>
    &emsp;0.5000000000&emsp;0.0000000000&emsp;0.0000000000&ensp;  20&ensp;  !&ensp;  M<br>
    &emsp;0.3333333333&emsp;0.3333333333&emsp;0.0000000000&ensp;  20&ensp;  !&ensp;  K<br>
    &emsp;0.0000000000&emsp;0.0000000000&emsp;0.0000000000&ensp;  20&ensp;  !&ensp;  G<br>
    &emsp;0.0000000000&emsp;0.0000000000&emsp;0.5000000000&ensp;  20&ensp;  !&ensp;  A<br>
    &emsp;0.5000000000&emsp;0.0000000000&emsp;0.5000000000&ensp;  20&ensp;  !&ensp;  L<br>
    &emsp;0.3333333333&emsp;0.3333333333&emsp;0.5000000000&ensp;  20&ensp;  !&ensp;  H<br>
    &emsp;0.0000000000&emsp;0.0000000000&emsp;0.5000000000&ensp;  20&ensp;  !&ensp;  A<br>
  
- inputなしでBZやkpathを出力するには、コード内の関数の`manual config`に結晶構造とkpathを直打ちし、<br>
ファイル名を引数に付けずに実行してください。<br>
<br>

## ssh先での表示
defaultでssh先かローカルかを判定し、ssh先ならばバックエンドをtkaggに変更して出力するようになっています。<br>
このため、<br>
&emsp;ssh -Y remote名<br>
でssh接続すればssh先でも使用できます。<br>
