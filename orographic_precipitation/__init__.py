from .orographic_precipitation import compute_orographic_precip

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("orographic_precipitation")
except PackageNotFoundError:  # noqa
    # package is not installed
    pass
