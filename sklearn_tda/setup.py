from distutils.core import setup

setup(
    name                           = "sklearn_tda",
    author                         = "Mathieu Carriere",
    author_email                   = "mathieu.carriere3@gmail.com",
    packages                       = ["sklearn_tda"],
    description                    = "A scikit-learn compatible package for doing Machine Learning and TDA",
    long_description_content_type  = "text/markdown",
    long_description               = open("README.md", "r").read(),
    url                            = "https://github.com/MathieuCarriere/sklearn_tda/",
    classifiers                    = ("Programming Language :: Python :: 3", "License :: OSI Approved :: MIT License", "Operating System :: OS Independent"),
)
