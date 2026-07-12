from setuptools import find_packages, setup


setup(
    name="coser",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
    ],
    entry_points={
        "console_scripts": [
            "coser=coser.cli:main",
        ],
    },
    author="shang",
    description="CoSEr coupling score calculator for paired traits.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.8",
)
