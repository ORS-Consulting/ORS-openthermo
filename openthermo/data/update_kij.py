import json
import os
from chemicals import *

"""
Script to add kij for n-alkanes / Water and set to 0.5 for all
"""


file = "pr.json"
name = "RAM-PR"
tables = {}
metadata = {}

file = os.path.join(os.path.abspath(os.path.dirname(__file__)), file)

f = open(file).read()
dat = json.loads(f)
data = dat["data"]

water_CAS = "7732-18-5"
HC_CASs = [
    search_chemical("methane").CASs,
    search_chemical("ethane").CASs,
    search_chemical("propane").CASs,
    search_chemical("i-butane").CASs,
    search_chemical("n-butane").CASs,
    search_chemical("i-pentane").CASs,
    search_chemical("n-hexane").CASs,
    search_chemical("heptane").CASs,
    search_chemical("n-octane").CASs,
    search_chemical("n-nonane").CASs,
    search_chemical("n-decane").CASs,
]
dict_template = {
    "Tmin": None,
    "P": None,
    "Pmin": None,
    "kij": 0,
    "T": None,
    "page": ["RAM"],
    "Pmax": None,
}

for CAS in HC_CASs:
    key = CAS + " " + water_CAS
    dict_template["kij"] = 0.5
    item = dict_template
    data[key] = item

import json

with open("pr_ram.json", "w") as fp:
    json.dump(dat, fp)
