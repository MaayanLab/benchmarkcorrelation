import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="benchmarkcorrelation",
    version="0.1.3",
    author="Alexander Lachmann",
    author_email="alexander.lachmann@mssm.edu",
    description="Benchmark correlation matrix for gene function prediction and general correctness.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/maayanlab/benchmarkcorrelation",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    package_data={
        "benchmarkcorrelation": ["data/*"]
    },
    include_package_data=True,
    install_requires=list(map(str.strip, open('requirements.txt', 'r').readlines())),
    python_requires='>=3.6',
)