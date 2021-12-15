from distutils.core import setup

setup(
    name='wave_sim',
    version='0.0.0.1',
    packages=["src.methods", "src.constants"],
    package_data={
         "": ["*"],
     },
    install_requires=[
        "click", "lhsmdu", "matplotlib", "numpy", "pyproj", "imageio", "tkinter"
    ],
)
