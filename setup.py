
from setuptools import setup, find_packages

setup(
    name="BZplot",
    version="2.0.2",
    author="Yudai Terao",
    author_email="terao.yudai.s4@dc.tohoku.ac.jp",
#    install_requires=[
#            'docopt',
#            'pathlib',
#            'typing',
#            'numpy',
#            'scipy',
#            'matplotlib',
#        ],
    entry_points={
        'console_scripts':[
            #'bz = BZplot.BZ:BZplot_main',
            'bz = BZplot.bzcell_exe:bz',
            #'cell = BZplot.BZ:cellplot_main',
            'cell = BZplot.bzcell_exe:cell',
        ],
    },
    packages=find_packages(),  # 使うモジュール一覧を指定する
)



