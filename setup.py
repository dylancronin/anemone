from setuptools import setup

setup(
    name="anemone",
    version="0.1.0",
    description="A standardized command-line pipeline integrating WGCNA and PLS-VIP regression",
    author="Dylan Cronin",
    url="https://github.com/dylancronin/anemone",
    scripts=["anemone"],
    python_requires=">=3.8",
    install_requires=[
        "pyyaml",
    ],
)
