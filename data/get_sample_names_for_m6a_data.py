import pandas as pd
import sys
import time

from pysradb import SRAweb 

SRP_list=["SRP007335", "SRP012098", "SRP012099" ,"SRP022152" ,"SRP026127" ,"SRP033229" ,"SRP039397", "SRP045349",
 "SRP059908" ,"SRP064809" ,"SRP067836" ,"SRP067910", "SRP076674", "SRP080113", "SRP082358" ,"SRP090254",
 "SRP090686", "SRP090687", "SRP094637", "SRP095597", "SRP096845", "SRP099081", "SRP099397", "SRP107200"]
db = SRAweb()
all_result=pd.DataFrame()
for srp in SRP_list:
        try:
            df = db.sra_metadata(srp)
            all_result=all_result.append(df)
        except:
            sys.stderr.write("Error with {}\n".format(srp))
            time.sleep(0.5)
        time.sleep(0.5)
all_result.to_csv("./data/all_srp.tsv", sep="\t", index=False)
