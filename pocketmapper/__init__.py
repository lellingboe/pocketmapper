"""
PocketMapper: compare binding pockets across structures via sequence or structural alignment.

Only `PocketMapper` and `__version__` are exported. The submodules are reachable as
`pocketmapper.lib.constants` etc. only as a side effect of `pocketmapper.pocketmapper` importing them, so
always import anything else explicitly -- `from pocketmapper.<module> import <name>`, or
`from pocketmapper.<package>.<module> import <name>` for a module in a subpackage.
"""

import logging

from pocketmapper.lib.constants import PACKAGE_LOGGER
from pocketmapper.pocketmapper import PocketMapper

# Library convention: records are dropped here unless the host application or the CLI adds a handler.
logging.getLogger(PACKAGE_LOGGER).addHandler(logging.NullHandler())

__all__ = ["PocketMapper", "__version__"]

__version__ = "0.2.6"
