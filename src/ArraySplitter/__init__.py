#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ArraySplitter")
except PackageNotFoundError:
    __version__ = "dev"

# __all__ = []
