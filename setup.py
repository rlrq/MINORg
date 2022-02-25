import setuptools

with open("README.md", 'r', encoding="utf-8") as fh:
    long_description = fh.read()

# Adapted from Biopython.
# We now define the MINORg version number in minorg/__init__.py
# Here we can't use "import minorg" then "minorg.__version__" as that would
# tell us the version of MINORg already installed (if any).
__version__ = "Undefined"
for line in open("minorg/__init__.py"):
    if line.startswith("__version__"):
        exec(line.strip())

setuptools.setup(
    name="minorg",
    version=__version__,
    author="Rachelle R.Q. Lee",
    author_email="rachelle.rq.lee@gmail.com",
    description="Generate minimum gRNA set for multiple non-reference genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rlrq/MINORg",
    project_urls={
        "Bug Tracker": "https://github.com/rlrq/MINORg/issues",
    },
    entry_points={
        "console_scripts": ["minorg=minorg.console:app"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"
    ],
    packages=[
        "minorg"
    ],
    include_package_data=True,
    install_requires=["biopython", "pybedtools", "pyfaidx", "regex", "typer"],
    python_requires=">=3.6",
)
