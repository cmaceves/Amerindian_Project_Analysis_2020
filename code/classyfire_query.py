import pandas as pd
import requests
import json

QUERY_ENDPOINT = "https://gnps-structure.ucsd.edu/classyfire?smiles="

def main():
    df = pd.read_csv("./classyfire_upload.csv")
    smiles = df["INCHI"].tolist()
   
    kingdom = []
    superclass = []
    class_list = []

    additional = []

    for smile in smiles:
        resp = requests.get(QUERY_ENDPOINT + str(smile))
        try:
            output = json.loads(resp.text)
        except:
            pass
        
        try:
            kingdom.append(output["kingdom"]["name"]) 
        except:
            kingdom.append("None")
        
        try:
            superclass.append(output["superclass"]["name"])
        except:
            superclass.append("None")
        try:
            class_list.append(output["class"]["name"])   
        except:
            class_list.append("None")

        
        try:
            additional_string = ""
            intermediates = output["intermediate_nodes"]
            for dictionary in intermediates:
                additional_string += dictionary["name"]
                additional_string += ","
            
            additional.append(additional_string)

        except:
            additional.append("None")


    df["kingdom"] = kingdom
    df["superclass"] = superclass
    df["class"] = class_list
    df["additional"] = additional
    df.to_csv("classyfire_output.csv")

            
if __name__ == "__main__":
    main()
