import setuptools

setuptools.setup(
    name="mol-ellipsize",
    version="0.0.1",
    author="Andrew Tarzia",
    author_email="andrew.tarzia@gmail.com",
    description="Ellipsoid fitting over conformer ensembles for calculating size.",
    url="https://github.com/andrewtarzia/mol-ellipsize",
    packages=setuptools.find_packages(),
    install_requires=(
        'scipy',
        'matplotlib',
        'networkx',
        'numpy',
    ),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
