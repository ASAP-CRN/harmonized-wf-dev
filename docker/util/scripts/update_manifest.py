#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd


def main(args):
    previous_data = pd.read_csv(args.previous_manifest, sep="\t")
    new_data = pd.read_csv(args.new_files_manifest, sep="\t")
    updated_data = (
        pd.concat([previous_data, new_data])
        .drop_duplicates("filename", keep="last")
        .sort_values(by="filename", axis=0)
    )
    updated_data.to_csv(args.updated_manifest, sep="\t", index=False)


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Merge file manifests; replace existing entries with updated information"
    )

    parser.add_argument(
        "-p",
        "--previous-manifest",
        type=str,
        required=True,
        help="Previous manifest to update",
    )
    parser.add_argument(
        "-n",
        "--new-files-manifest",
        type=str,
        required=True,
        help="Set of new or updated files",
    )
    parser.add_argument(
        "-u",
        "--updated-manifest",
        type=str,
        required=True,
        help="Output file to write updated manifest to",
    )

    args = parser.parse_args()
    main(args)
