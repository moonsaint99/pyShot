from setuptools import setup, find_packages

setup(
    name="pyshot",
    version="0.0.1",
    packages=find_packages(),
    install_requires=[
        "matplotlib",
        "obspy",
        "scipy",
        "numpy",
        "pandas",
        "pyqt5"
        # Add other dependencies here
    ],
    author="Rishi Chandra",
    author_email="chandra@arizona.edu",
    description="A tool for picking seismic events in a seismic signal",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/moonsaint99/pyShot",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
