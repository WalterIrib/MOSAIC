from setuptools import setup, find_packages

setup(
    name="mosaic-pipeline",
    version="1.1.0",
    description="Multi-layer Omni-transcriptomic Sequencing Analysis and Integrated Characterization",
    author="Walter Arancio",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "click>=8.0",
        "PyYAML>=6.0",
        "requests>=2.28",
    ],
    entry_points={
        "console_scripts": [
            "mosaic=mosaic.cli:cli",
            "mosaic-refs=mosaic.cli:refs_cli",
            "mosaic-run=mosaic.cli:run_cli",
        ],
    },
)
