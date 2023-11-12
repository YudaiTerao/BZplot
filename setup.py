
from setuptools import setup

setup(
    name="BZplot",
    version="0.2.0",
    author="Yudai Terao",
    entry_points={
        'console_scripts':[
            'bz = BZplot.BZ:main',
        ],
    },
)



