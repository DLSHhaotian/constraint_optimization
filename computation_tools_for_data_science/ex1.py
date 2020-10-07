# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 17:07:58 2020

@author: DLSH
"""

import pandas as pd 
import numpy as np
import re
import ast
import math
#import simplejson as js #not used use the eval
url = 'http://courses.compute.dtu.dk/02807/2020/projects/project1/movies_metadata.csv'

def load_movies_data(url):
    df = pd.read_csv(url, encoding='utf-8')
    
    #use the in-build string function .str()
    #del var (delete)
    #use the in-build string function .to_datetime from str to TimeStamped
    df=df[(df.release_date.str.len())>5]#already delete the bad data
    df.release_date = pd.to_datetime(df.release_date, format = '%Y-%m-%d')
    #eval(var) from string to dict
    #type(var) return the type of var
    #apply(func) = loop(func); but return a new one  
    str2dict_rows=df.belongs_to_collection.apply(type)==str
    df.update(df[str2dict_rows].belongs_to_collection.apply(ast.literal_eval))
    
    df.update(df.genres.apply(ast.literal_eval))
    df.update(df.production_companies.apply(ast.literal_eval))
    df.update(df.production_countries.apply(ast.literal_eval))
    return df
    
import time
start = time.time()
df=load_movies_data(url)
end = time.time()
print (str(end-start))

df=df[df.adult.apply(lambda x:True if  x=='False' else False)]#drop if adult is true

