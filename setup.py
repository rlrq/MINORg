import setuptools

with open("README.md", 'r', encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MINORg_rlrq",
    version="2.2",
    author="Rachelle R.Q. Lee",
    author_email="rachelle2nd@yahoo.com.sg",
    description="Generate minimum gRNA set for multiple non-reference genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rlrq/MINORg",
    project_urls={
        "Bug Tracker": "https://github.com/rlrq/MINORg/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"
    ],
    package_dir={"": "minorg"},
    packages=setuptools.find_packages(where="minorg")
    python_requires=">=3.6",
)
