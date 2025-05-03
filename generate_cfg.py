import yaml
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--config", help="Base configfile.")
    ap.add_argument()
    args = ap.parse_args()

    pass


if __name__ == "__main__":
    raise SystemExit(main())
