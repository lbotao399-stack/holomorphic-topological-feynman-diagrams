"""Replay the complete supported-source polynomial certificate and controls."""
import argparse
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deep", action="store_true", help="regenerate all complete source samples and controls")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    if args.deep:
        with TemporaryDirectory(prefix="replay_", dir=HERE) as temporary:
            subprocess.run([sys.executable, str(HERE / "certify.py"), "--jobs", str(args.jobs),
                            "--checkpoint", temporary], check=True)
        subprocess.run([sys.executable, str(HERE / "independent_controls.py"), "--recompute"], check=True)
    subprocess.run([sys.executable, str(HERE / "verify.py")], check=True)
    subprocess.run([sys.executable, str(HERE / "render.py")], check=True)


if __name__ == "__main__":
    main()
