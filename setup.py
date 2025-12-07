from setuptools import setup, find_packages

setup(
    name="bio_scripts",
    version="0.1.0",
    description="Bioinformatics analysis scripts and utilities",
    package_dir={"": "src/python"},
    packages=find_packages(where="src/python"),
    python_requires=">=3.8",
    install_requires=[
        "pandas",
        "numpy",
        # Add other dependencies here
    ],
)
