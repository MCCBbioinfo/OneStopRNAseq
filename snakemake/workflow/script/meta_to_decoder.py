import pandas as pd
import sys

def read_table(fname='meta/contrast.de.xlsx'):
    try:
        if fname.endswith(".txt"):
            df = pd.read_table(fname, keep_default_na=False, na_values=[])
        elif fname.endswith(".csv"):
            df = pd.read_csv(fname, keep_default_na=False, na_values=[])
        elif fname.endswith(".xlsx"):
            df = pd.read_excel(fname, engine='openpyxl', keep_default_na=False, na_values=[])
        else:
            sys.exit("fname not xlsx nor txt")
        df = df.replace("NA", "NA_")
    except:
        sys.exit(">>> Fatal Error: " + fname + " format error, it can't be read correctly." +
                 "\nPlease check if you saved txt file as xlsx or vice versa\n\n")
    return (df)

fname = sys.argv[1]
df = read_table(fname)
df.columns = ["sample.ID", "group.ID", "batch.ID"]
df.to_csv("meta/decoder.txt", sep="\t", index=False)
