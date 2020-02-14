import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name="orographic_precipitation",
     version='0.1a',
     description="Linear Theory of Orographic Precipitation",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/rlange2/orographic-precipitation",
     download_url="https://github.com/rlange2/orographic-precipitation/archive/0.1a.tar.gz",
     license = "MIT",
     maintainer="Raphael Lange",
     maintainer_email="raphael.lange@gfz-potsdam.de",
     packages=setuptools.find_packages(),
     python_requires=">=3.6",
     install_requires=["matplotlib", "numpy"]
 )
