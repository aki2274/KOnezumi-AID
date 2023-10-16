import re
from typing import Tuple
import numpy as np
import pandas as pd

def set_presetnum(name:str,gene_df):
    data = gene_df[gene_df['name']==name]
    cdsStart =int(data['cdsStart'])
    start = data['exonStarts']
    end = data['exonEnds']
    start_list = list(start.split(','))
    end_list = list(end.split(','))
    return data,cdsStart,start_list,end_list

