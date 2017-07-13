from __future__ import print_function, absolute_import
import sys
from hic2cool.__main__ import main

warning = """The `run_hic2cool.py` script is deprecated and will be removed. 
Instead, please use `python -m hic2cool` or just `hic2cool` if the hic2cool 
package is installed in your environment.
"""

print(warning, file=sys.stderr)
main()
