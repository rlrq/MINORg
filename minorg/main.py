#!/usr/bin/env python3

"""
Calls MINORg app.
"""

import os
import sys

## just cheat it a little so the minorg.X imports work without this actually being installed
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

# from minorg.console import app
from minorg.console import app

app()
