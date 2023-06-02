from setuptools import setup
#from setuptools.extension import Extension
#from Cython.Build import cythonize

#extensions = [
#    Extension("st_rle.rle", ["st_rle/rle.py"])
#]

setup(
    name="st_rle",
    version="0.1.0",
    author="Scott Tyler",
    author_email="scottyler89@gmail.com",
    url="https://github.com/scottyler89/st_rle",
    description="RLE normalization",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=["st_rle"],
    install_requires=[
        "numpy",
        "scipy",
        "joblib",
        "cython",
    ],
    #ext_modules=cythonize(extensions),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)

