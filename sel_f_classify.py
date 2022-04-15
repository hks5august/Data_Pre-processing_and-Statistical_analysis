# Feature Extraction with Univariate Statistical Tests (ANOVA-F values for classification)

# usage:
# python sel_f_classify.py <input file WITHOUT HEADER having comma-separated features, last column label 1,-1 or 1,0> <output file> <No. of features to be selected> 

import sys
import pandas
import numpy as np
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif

infile = sys.argv[1]
outfile = sys.argv[2]
f = open(outfile,'a')
featnum = int(sys.argv[3])


# load data
#url = "https://archive.ics.uci.edu/ml/machine-learning-databases/pima-indians-diabetes/pima-indians-diabetes.data"
#names = ['preg', 'plas', 'pres', 'skin', 'test', 'mass', 'pedi', 'age', 'class']
df = pandas.read_csv(infile)
Y  = df['flag']
X = df.ix[:, df.columns != "flag"]
# feature extraction
test = SelectKBest(score_func=f_classif, k=featnum)
fit = test.fit(X, Y)
#print(fit.get_support())
a=np.where(fit.get_support() == True)[0]
#print(a)
head= (list(df))
for n in a:
	#print(type(n))
	f.write(head[n])
	f.write(",")
