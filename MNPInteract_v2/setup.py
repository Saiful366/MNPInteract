from setuptools import setup, find_packages

setup(
    name="MNPInteract",
    version="1.1.0",
    author="Saiful Islam",
    description=(
        "Complete post-AlphaPulldown pipeline for identifying high-confidence "
        "PDLP5-interacting proteins using score filtering, PFAM/GO annotation, "
        "DeepTMHMM topology prediction, atom-level interface detection, "
        "and high-confidence interactor identification."
    ),
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "pandas",
        "biopython",
        "requests",
        "pybiomart",
    ],
    entry_points={
        "console_scripts": [
            "MNPInteract=MNPInteract.MNPInteract:main",
        ]
    },
)
