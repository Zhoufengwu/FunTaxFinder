"""FunTaxFinder package."""

from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version("funtaxfinder")
except PackageNotFoundError:
    __version__ = "0+unknown"


__all__ = ["__version__"]
