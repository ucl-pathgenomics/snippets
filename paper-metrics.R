# req R & Python

library(reticulate)
library(tidyverse)

# Seeing your enviroments
conda_list()


#Using it
conda_list()[[1]][2] %>%
  use_condaenv(required = TRUE)

repl_python()

import requests
from bs4 import BeautifulSoup as bs # HTML parser
import os # Useful for file system operations
import datetime

date = datetime.datetime.now()
date = date.strftime("%Y-%m-%d")

URL = 'https://www.biorxiv.org/content/10.1101/2020.05.15.097907v1.article-metrics'
header ={'User-Agent': 'Mozilla/5.0 (Windows NT x.y; Win64; x64; rv:10.0) Gecko/20100101 Firefox/10.0 '}
r = requests.get(URL, headers=header)

soup = bs(r.content, 'html.parser')

results = soup.find("tr",{"class":"odd"})

headers = ['date', 'ignore', 'abstract', 'full', 'pdf']
data = [date]

for row in results:
  data.append(row.string)


outfile = open("/home/oscar/regularjobs/paper-stats/data.tab", "a")
outfile.write("\n") \
outfile.write("\t".join(data))

exit
