import setuptools


def open_requirements(path):
    with open(path) as f:
        requires = [r.split("/")[-1] if r.startswith("git+") else r for r in f.read().splitlines()]
    return requires


# Get requirements from requirements.txt, stripping the version tags
requires = open_requirements("requirements.txt")

setuptools.setup(
    name="cylinterp",
    version="0.0.1",
    packages=setuptools.find_packages(),
    package_data={'cylinterp': ['data_files/*']},
    install_requires=requires,
    description="Cylindrical interpolation package for drifting electrons",
    author="Jianyang Qi",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    zip_safe=False,
)