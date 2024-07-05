
from setuptools import setup, find_packages

setup(
    name="BZplot",
    version="1.2.2",
    author="Yudai Terao",
    entry_points={
        'console_scripts':[
            'bz = BZplot.BZ:BZplot_main',
            'cell = BZplot.BZ:cellplot_main',
        ],
    },
    packages=find_packages(),  # 使うモジュール一覧を指定する
)



