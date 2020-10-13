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
#apply(axis=0,column,axis=1,raw)
#define a function which can judge the NAN and return the (float)1 or 0 for next sum step 
def rating_count(x):
    return x.apply(lambda x:0.0 if np.isnan(x) else 1.0)
#Traverse all return the 1 or 0
rating_count_float=user_ratings.apply(rating_count,axis=0)
#sum rating count for each column(movie)
rating_count_sum=rating_count_float.apply(sum,axis=0)
#convert to a dataframe and rename the column
rating_count_sum=pd.DataFrame(rating_count_sum).rename(columns={0:'rating_counts'})
#merge with the 'title'
toy_corr_new=pd.merge(toy_corr, rating_count_sum, on='title',how='left')
#delete Toy Story (1995) itself
toy_corr_new=toy_corr_new[~toy_corr_new.title.str.contains(r'Toy Story \(1995\)')].reset_index()
#define a function to find wanted rating counts' movies
def rating_count_threshold(x):
    return True if x>100 else False
toy_corr_new[toy_corr_new.rating_counts.apply(rating_count_threshold)].sort_values(by=['Correlation'],ascending=False).head(5)
