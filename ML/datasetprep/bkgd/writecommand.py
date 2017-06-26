import numpy as np
import pickle
import os

ids = pickle.load(open("out.pickle", "rb"))

countids = 0
countfiles = 1

output = open("commandbkgd1.sql", "w+")
output.write("select snobjid,ra,dec,mjd from WSDIFF.SNOBS where snobjid in (")
for aid in ids:
    if countids <= 950:
        countids = countids + 1
        if countids == 951:
            output.write(str(aid))
        else:
            output.write(str(aid) + ', ')
    else:
        output.write('); > mjdbkgd' + str(countfiles) + ".csv")
        countids = 1
        countfiles= countfiles + 1
        output = open("commandbkgd"+str(countfiles) + ".sql", "w+")
        output.write("\n")
        output.write("select snobjid,ra,dec,mjd from WSDIFF.SNOBS where snobjid in (")
        output.write(str(aid)+', ')
  
output.write('); > mjdbkgd' + str(countfiles) + ".csv")
