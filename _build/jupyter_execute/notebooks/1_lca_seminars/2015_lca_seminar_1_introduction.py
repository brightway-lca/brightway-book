#!/usr/bin/env python
# coding: utf-8

# # Class 1: Introduction and setting up a simple LCA calculation
# 
# Note: this notebook is based on [Brightway2 tutorial 4 - meta-analysis](http://nbviewer.ipython.org/urls/bitbucket.org/cmutel/brightway2/raw/default/notebooks/Tutorial%204%20-%20Meta-analysis.ipynb).
# 
# # Setup
# 
# Import the necessary libraries

# In[1]:


from brightway2 import *
import numpy as np
import pyprind
from scipy import stats
import random


# Set a new project for this class session

# In[2]:


projects.current = "Class 1"


# Import the basic biosphere and LCIA methods (requires internet connection)

# In[3]:


bw2setup()


# # Import ecoinvent 2.2
# 
# We are using version 2.2 to make the calculations a bit quicker, as this is intended to be run during the class.
# 
# First, extract the data from the XML files.

# In[4]:


ei = SingleOutputEcospold1Importer("/Users/cmutel/Documents/LCA Documents/Ecoinvent/2.2/processes", "ecoinvent 2.2")


# On windows, you will need to escape your backslashes by repeating them twice, e.g.

# In[ ]:


ei = SingleOutputEcospold1Importer("c:\\Users\\cmutel\\Process_infra_roh", "ecoinvent 2.2")


# Next, we normalize some values, and link the different datasets to each other and the basic biosphere

# In[5]:


ei.apply_strategies()


# Check to make sure everything is linked

# In[6]:


ei.statistics()


# Finally, write the database

# In[7]:


ei.write_database()


# In[9]:


list(databases)


# In[10]:


list(methods)[:10]


# # Do some calculations
# 
# We select the new database we just created, and then use a bit of trickiness to get random ordering. All datasets have the ``type`` ``"process"``, so this is equivalent to random ordering.

# In[11]:


db = Database("ecoinvent 2.2")
db.order_by = "type"


# Many of the built-in LCIA methods are too narrowly focused - we pre-select a set of candidates to use for our meta-calculations.

# In[12]:


CANDIDATES = sorted([
 (u'CML 2001', u'acidification potential', u'average European'),
 (u'CML 2001', u'climate change', u'GWP 100a'),
 (u'CML 2001', u'eutrophication potential', u'average European'),
 (u'CML 2001', u'freshwater aquatic ecotoxicity', u'FAETP 100a'),
 (u'CML 2001', u'human toxicity', u'HTP 100a'),
 (u'CML 2001', u'land use', u'competition'),
 (u'CML 2001', u'marine aquatic ecotoxicity', u'MAETP infinite'),
 (u'CML 2001', u'resources', u'depletion of abiotic resources'),
 (u'CML 2001', u'stratospheric ozone depletion', u'ODP 25a'),
 (u'EDIP2003', u'ecotoxicity', u'in sewage treatment plants'),
 (u'EDIP2003', u'eutrophication', u'terrestrial eutrophication'),
 (u'EDIP2003', u'renewable resources', u'wood'),
 (u'EDIP2003', u'stratospheric ozone depletion', u'ODP total'),
 (u'EPS 2000', u'total', u'abiotic stock resources'),
 (u'EPS 2000', u'total', u'emissions into soil'),
 (u'EPS 2000', u'total', u'emissions into water'),
 (u'EPS 2000', u'total', u'land occupation'),
 (u'IMPACT 2002+ (Endpoint)', u'ecosystem quality', u'land occupation'),
 (u'IMPACT 2002+ (Endpoint)', u'human health', u'ozone layer depletion'),
 (u'IMPACT 2002+ (Endpoint)', u'resources', u'mineral extraction'),
 (u'IMPACT 2002+ (Endpoint)', u'resources', u'non-renewable energy'),
 (u'IMPACT 2002+ (Midpoint)', u'ecosystem quality', u'aquatic acidification'),
 (u'IPCC 2001', u'climate change', u'GWP 100a'),
 (u'ReCiPe Endpoint (H,A)',
  u'ecosystem quality',
  u'agricultural land occupation'),
 (u'ReCiPe Endpoint (H,A)',
  u'ecosystem quality',
  u'freshwater eutrophication'),
 (u'ReCiPe Endpoint (H,A)',
  u'ecosystem quality',
  u'natural land transformation'),
 (u'ReCiPe Endpoint (H,A)',
  u'ecosystem quality',
  u'terrestrial acidification'),
 (u'ReCiPe Endpoint (H,A)', u'ecosystem quality', u'urban land occupation'),
 (u'ReCiPe Endpoint (H,A)', u'human health', u'particulate matter formation'),
 (u'ReCiPe Endpoint (H,A)', u'resources', u'fossil depletion'),
 (u'TRACI', u'environmental impact', u'acidification'),
 (u'TRACI', u'environmental impact', u'eutrophication'),
 (u'TRACI', u'environmental impact', u'global warming'),
 (u'TRACI', u'environmental impact', u'ozone depletion'),
 (u'TRACI', u'human health', u'respiratory effects, average'),
 (u'eco-indicator 99, (H,A)',
  u'ecosystem quality',
  u'acidification & eutrophication'),
 (u'eco-indicator 99, (H,A)', u'ecosystem quality', u'ecotoxicity'),
 (u'eco-indicator 99, (H,A)', u'ecosystem quality', u'land occupation'),
 (u'eco-indicator 99, (H,A)', u'human health', u'carcinogenics'),
 (u'eco-indicator 99, (H,A)', u'human health', u'climate change'),
 (u'eco-indicator 99, (H,A)', u'human health', u'ozone layer depletion'),
 (u'eco-indicator 99, (H,A)', u'resources', u'fossil fuels'),
 (u'eco-indicator 99, (H,A)', u'resources', u'mineral extraction'),
 (u'ecological footprint', u'total', u'CO2'),
 (u'ecological footprint', u'total', u'land occupation'),
 (u'ecological footprint', u'total', u'nuclear'),
 (u'ecological scarcity 2006', u'total', u'deposited waste'),
 (u'ecological scarcity 2006', u'total', u'emission into groundwater'),
 (u'ecological scarcity 2006', u'total', u'energy resources'),
 (u'ecological scarcity 2006', u'total', u'natural resources'),
 (u'ecosystem damage potential', u'total', u'linear, land occupation'),
 (u'ecosystem damage potential', u'total', u'linear, land transformation'),
])

assert all(x in methods for x in CANDIDATES)

print("There are %s methods to test" % len(CANDIDATES))


# Choose ten LCIA methods and 500 datasets from ecoinvent 2.2 at random

# In[13]:


chosen_methods = random.sample(CANDIDATES, 10)
chosen_processes = []
for index, obj in enumerate(db):
    if index >= 500:
        break
    else:
        chosen_processes.append(obj)


# Set up the LCA object, optimized to do many calculations.
# 
# See [making LCA calculations faster](http://chris.mutel.org/fast-dont-lie.html) blog post for more details on factorization.

# In[14]:


lca = LCA({chosen_processes[0]: 1}, method=chosen_methods[0])
lca.lci(factorize=True)
lca.lcia()


# Create an array to store our LCA results - processes on rows, methods on columns

# In[15]:


results = np.zeros((500, 10))


# Do 5000 LCA calculations in a single thread. Store the results in ``results``.

# In[16]:


bar = pyprind.ProgBar(5000, monitor=True)

for col, method in enumerate(chosen_methods):
    lca.method = method
    lca.load_lcia_data()
    for row, process in enumerate(chosen_processes):
        lca.redo_lcia({process: 1})
        results[row, col] = lca.score
        bar.update()

print(bar)


# We only care about processes which have non-zero LCA scores - there are a few processes in ecoinvent 2.2 which we want to filter automatically (if they are selected in our earlier random sample).

# In[17]:


mask = (results.sum(axis=1) != 0)
print("Ignoring {} processes".format((~mask).sum()))


# Calculate the rank-order correlation for all processes

# In[18]:


def create_correlation_matrix(scores_array):
    num_methods = scores_array.shape[1]
    correlations = np.zeros((num_methods, num_methods))

    for row in range(num_methods):
        for col in range(num_methods):
            if col <= row:
                continue                               # Only need to compute correlation once
            dataset_1 = scores_array[:, row]
            dataset_2 = scores_array[:, col]
            mask = (dataset_1 != 0) * (dataset_2 != 0) # Ignore activities that have zero score
            correlations[row, col] = stats.kendalltau( # Get tau value, drop p-statistic
                dataset_1[mask], 
                dataset_2[mask]
            )[0]  

    return correlations


# In[19]:


correlation_matrix = create_correlation_matrix(results[mask, :])


# Visualize the results

# In[20]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[21]:


import matplotlib.pyplot as plt

fig = plt.gcf()
fig.set_size_inches(12, 12)

masked_correlation = np.ma.array(correlation_matrix, mask=correlation_matrix == 0).T
plt.pcolor(masked_correlation, cmap=plt.cm.cubehelix_r)
plt.colorbar(label=r"Kendall $\tau$ rank-order correlation coefficient")
plt.ylim(None, 10)
plt.xlim(None, 10)
plt.axis('off')

