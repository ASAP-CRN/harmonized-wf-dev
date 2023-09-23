import pandas as pd
from pathlib import Path

## load tables
metadata_path = Path.cwd()/ "metadata"

SUBJECT = pd.read_csv(f"{metadata_path}/SUBJECT.tsv", delimiter="\t")
SAMPLE = pd.read_csv(f"{metadata_path}/SAMPLE.tsv",delimiter="\t")
CLINPATH = pd.read_csv(f"{metadata_path}/CLINPATH.csv",delimiter=",")

STUDY = pd.read_csv(f"{metadata_path}/STUDY.tsv",delimiter="\t")
PROTOCOL = pd.read_csv(f"{metadata_path}/PROTOCOL.tsv",delimiter="\t")

# fix STUDY 
tmp = pd.DataFrame()
tmp = STUDY[["Unnamed: 1","Unnamed: 0"]].transpose().reset_index().drop(columns=["index"])
tmp.columns = tmp.iloc[0]
STUDY = tmp.drop([0])

# write clean tables to metadata/clean/
STUDY.to_csv(metadata_path / "clean/STUDY.csv")
PROTOCOL.to_csv(metadata_path / "clean/PROTOCOL.csv")
CLINPATH.to_csv(metadata_path / "clean/CLINPATH.csv")
SAMPLE.to_csv(metadata_path / "clean/SAMPLE.csv")
SUBJECT.to_csv(metadata_path / "clean/SUBJECT.csv")
