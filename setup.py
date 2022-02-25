import setuptools

with open("README.md", 'r', encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="minorg",
    version="0.2.1rc",
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
