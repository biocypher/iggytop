import pandas as pd


class IEDBAdapter:

    def __init__(self):
        table = pd.read_csv("data/iedb_test.csv")
        column_names = table.iloc[0].to_dict()
        column_rename = {}
        for key, value in column_names.items():
            column_rename[key] = key.split(".")[0] + " " + value

        table = table.rename(columns=column_rename)
        table = table.drop(0)

        self.epitopes = table[[
            "Epitope Name", 
            "Epitope Source Molecule", 
            "Epitope Source Organism"
        ]].drop_duplicates()
