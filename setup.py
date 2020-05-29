import setuptools

setuptools.setup(
    name="cellscape",
    version="0.0.0",
    packages=setuptools.find_namespace_packages(),
    entry_points={'console_scripts':['cellscape = cellscape.__main__:main']}
)
