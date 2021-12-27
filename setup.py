from distutils.core import setup

setup(
    name='wave_sim',
    version='0.0.0.3',
    packages=["src.schrodinger", "src.analysis", "src.gui"],
    package_data={
         "": ["*"],
     },
    install_requires=[
        "click", "pandas", "matplotlib", "numpy", "pyproj", "imageio", "tkinter"
    ],
)
