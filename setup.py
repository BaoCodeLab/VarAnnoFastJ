from setuptools import setup, find_packages

setup(
    name="compGWAS",
    version='1.0',
    description="A python tool for variant annotation.",
    classifiers=['License :: GPL-3.0','Operating System :: POSIX :: Linux'],
    author="Yun-Juan Bao",
    author_email="yjbao@hubu.edu.cn"
    url="https://github.com/BaoCodeLab/VarAnnoFastJ",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=["numpy >= 1.17.3","pandas >= 1.1.4","gffutils >= 0.10.1"]



