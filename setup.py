from setuptools import setup, find_packages

setup(
    name="python_script",
    version="0.1.0",
    description="Analysis scripts and utilities by dazheng",
    package_dir={"": "src/python"},
    packages=find_packages(where="src/python"),
    python_requires=">=3.12",
    install_requires=[
        "cookiecutter",
        "fastcluster",
        "matplotlib",
        "numpy",
        "pandas",
        "PyYAML",
        "scikit-learn",
        "scipy",
        "seaborn",
        "upsetplot",
        "statsmodels",
    ],
)
