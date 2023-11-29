
from setuptools import setup

setup(
    name="BZplot",
    version="0.3.1",
    author="Yudai Terao",
    entry_points={
        'console_scripts':[
            'bz = BZplot.BZ:main',
        ],
    },
)



