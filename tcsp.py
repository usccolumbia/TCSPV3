"""Local TCSPV3 entry point."""
import sys
from pathlib import Path
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
from tcspv3.__main__ import main
if __name__ == '__main__':
    main()
