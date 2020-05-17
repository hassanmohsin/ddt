from ddt.utility import SmilesToImage, SDFToImage


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert SMILES to PNG")
    input = parser.add_mutually_exclusive_group(required=True)
    input.add_argument(
        "--smiles", action="store", dest="smiles", required=False, help="SMILES string",
    )
    input.add_argument(
        "--sdf", action="store", dest="sdf", help="SDF file",
    )
    parser.add_argument(
        "--out", action="store", dest="out", required=True, help="PNG file location"
    )
    args = parser.parse_args()

    if args.smiles:
        stoi = SmilesToImage(args.smiles)
        stoi.toPNG(output=args.out)
        # stoi.toPNG(output=args.out, get_binary=True)
    if args.sdf:
        stoi = SDFToImage(args.sdf)
        stoi.toPNG(output=args.out, max=8)
        # stoi.toPNG(output=args.out, get_binary=True)

