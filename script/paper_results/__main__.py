"""Explicit recompute / render commands."""
import argparse
from pathlib import Path

def main():
    parser=argparse.ArgumentParser(description="Paper03 geometric sectors and figure supplement.")
    sub=parser.add_subparsers(dest="command",required=True)
    a=sub.add_parser("recompute");a.add_argument("--output-dir",required=True)
    b=sub.add_parser("render");b.add_argument("--output-dir",required=True)
    b.add_argument("--data-dir",default=str(Path(__file__).resolve().parents[2]/"data/revised"))
    args=parser.parse_args()
    if args.command=="recompute":
        from .data import recompute
        print(recompute(args.output_dir))
    else:
        from .figures import render
        for path in render(args.data_dir,args.output_dir):print(path)

if __name__=="__main__":main()
