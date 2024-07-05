
from setuptools import setup

setup(
    name="BZplot",
    version="1.2.0",
    author="Yudai Terao",
    entry_points={
        'console_scripts':[
            'bz = BZplot.BZ:BZplot_main',
            'cell = BZplot.BZ:cellplot_main',
        ],
    },
)



