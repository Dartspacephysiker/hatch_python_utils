# 2019/10/02
import re

# myRE = re.compile("NCAR_([0-9]{4})-([0-9]{2})-([0-9]{2})_([a-zA-Z0-9]+)_([a-zA-Z0-9.]+)_([a-zA-Z0-9]+)")
# myRE = re.compile("([0-9]{4})-([0-9]{2})-([0-9]{2})_([a-zA-Z0-9]+)__?([a-zA-Z0-9]+)_([a-zA-Z0-9.]+).pkl")
# iHere = 0
# for i,samedate in enumerate(samedates):
#     matcha = myRE.search(samedate)
#     if matcha is not None:
#         yr,mo,day,radnavn,style,integtid = matcha.groups()
#         dateDict[iHere] = dict(yr=yr,day=day,mo=mo,radnavn=radnavn,style=style,integtid=integtid)
#         iHere += 1
