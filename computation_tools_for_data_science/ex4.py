# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 13:55:56 2020

@author: DLSH
"""

import pandas as pd 
import numpy as np
import re
import ast
import math

url1 = 'http://courses.compute.dtu.dk/02807/2020/projects/project1/ratings.csv'
url2 = 'http://courses.compute.dtu.dk/02807/2020/projects/project1/movies.csv'

df1 = pd.read_csv(url1, encoding='utf-8')
df2 = pd.read_csv(url2, encoding='utf-8')
#similiar to the join in sql just connect two dataframe
movie_data=pd.merge(df1, df2, on='movieId',how='left')
#pick the used info into a new dataframe
movie_data_rating=pd.DataFrame(movie_data,columns=(['userId','rating','title']))
#using pivot_table to set specific column or index
user_ratings = movie_data_rating.pivot_table(columns=['title'], values='rating', index='userId')
#pick the toystory info
toystory_ratings = user_ratings['Toy Story (1995)']
#cal the correlation
toy_corr=pd.DataFrame(user_ratings.corrwith(toystory_ratings,method='pearson'))
toy_corr=toy_corr.reset_index()
toy_corr=toy_corr.rename(columns={0:'Correlation'})
toy_corr.sort_values(by=['Correlation'],ascending=False).head(5)

def rating_count(x):
