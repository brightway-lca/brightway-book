#!/usr/bin/env python
# coding: utf-8

# # Class 4 - Hybrid LCA
# 
# In this class, we will learn about supply use tables, and input output tables. We will also do a toy hybrid LCA.
# 
# Before getting started, make sure you have upgrade the Brightway2 packages. You should have at least the following:

# In[1]:


import bw2data, bw2calc, bw2io
print("BW2 data:", bw2data.__version__)
print("BW2 calc:", bw2calc.__version__)
print("BW2 io:", bw2io.__version__)


# Now import the necessary libraries:

# In[2]:


from brightway2 import *
from bw2io.importers.exiobase import Exiobase22Importer
import numpy as np
import os
import pyprind


# Create a new project for this class:

# In[3]:


if 'Class 4' not in projects:
    projects.current = "Class 1"
    projects.copy_project("Class 4")

projects.current = "Class 4"


# We will need the latest version of the data migrations to match EXIOBASE biosphere flows to ecoinvent biosphere flows:

# In[4]:


create_core_migrations()


# In[5]:


ERROR_MSG = """Missing a data migration needed for this class. 

Please make sure you hvae the latest Brightway2 libraries, and reset the notebook."""
assert 'exiobase-biosphere' in migrations, ERROR_MSG


# # Import EXIOBASE 2.2
# 
# Now we need to download the industry by industry table from version 2.2 of exiobase. You can get it from the following link. Note that you will have to register an account if this is the first time you use this database: http://www.exiobase.eu/index.php/data-download/exiobase2-year-2007-full-data-set/78-mriot-ixi-fpa-coefficient-version2-2-2/file
# 
# Extract the downloaded file, and adjust the following. Windows users might need something like:
# 
#     fp = "C:\\Users\\<your name>\\Downloads\\mrIOT_IxI_fpa_coefficient_version2.2.2"

# In[6]:


fp = "/Users/cmutel/Downloads/mrIOT_IxI_fpa_coefficient_version2.2.2"

assert os.path.exists(fp), "Please adjust your filepath, the provided one doesn't work"


# We can now import the exiobase database. This will take a while, so go ahead and get started.
# 
# Why is this so slow compared to ecoinvent, for example? The answer lies in the density of the technosphere matrix. Exiobase, and IO tables in general, use comprehensive data from surveys and national customs, so they will get data on things that normal people would never even think of. For example, how much rice from Thailand is required to produce one euro of steel in Germany?
# 
# In other words, the technosphere matrix is very dense. Ecoinvent is stored as a [sparse matrix](http://docs.scipy.org/doc/scipy/reference/sparse.html), where data is only provided in about 1.5% of all possible locations - every other value is zero, and these zeros are not stored, only implied. However, the IO table has a fill rate of about 50%, meaning that we store every value in the matrix. The technosphere in ecoinvent 2.2 is about 4000 by 4000, but we only need to store about 40.000 numbers. The technosphere matrix is exiobase is about 8000 by 8000, but we store around 35.000.000 numbers.
# 
# We use a special backend for IO databases, as our standard storage mechanisms simply fall apart with such large data sets. You can see this [backend here](https://bitbucket.org/cmutel/brightway2-data/src/tip/bw2data/backends/iotable/__init__.py?at=2.0&fileviewer=file-view-default).

# In[7]:


ex = Exiobase22Importer(fp)
ex.apply_strategies()
ex.write_database()


# Free up some memory

# In[8]:


ex = None


# # LCA calculations
# 
# We can now do an LCA. We first do this the standard way:

# In[ ]:


gwp = ('IPCC 2013', 'climate change', 'GWP 100a')
lca = LCA({Database("EXIOBASE 2.2").random(): 1}, method=gwp)
lca.lci()
lca.lcia()


# Our technosphere matrix is sparse:

# In[13]:


lca.technosphere_matrix


# And it takes a while to solve (versus less than one second for ecoinvent 2.2):

# In[14]:


get_ipython().run_line_magic('timeit', 'lca.solve_linear_system()')


# Free up some memory by forgetting about the `lca` object.

# In[15]:


lca = None


# However, we have a special LCA class that only does [dense technosphere matrices](https://bitbucket.org/cmutel/brightway2-calc/src/tip/bw2calc/dense_lca.py?at=default&fileviewer=file-view-default). If we use it, we will get better performance, because the linear solver assumes dense instead of sparse matrices:

# In[16]:


dlca = DenseLCA({Database("EXIOBASE 2.2").random(): 1}, method=gwp)
dlca.lci()


# The technosphere is, as you would expect, now a dense matrix

# In[17]:


type(dlca.technosphere_matrix)


# The nupy dense solver of linear system is faster than the SciPy/UMFPACK sparse solver, as our matrix actually is quite dense. The performance should be much better:

# In[18]:


get_ipython().run_line_magic('timeit', 'dlca.solve_linear_system()')


# Free up some more memory by forgetting about the `tech_params` array.

# In[19]:


print(dlca.tech_params.shape)
dlca.tech_params = None


# # Create aggregated processes
# 
# We can now create aggregated (so-called "system") processes for each activity in Exiobase. These aggregated proceses can be used in our normal sparse LCAs, but are terminated, i.e. we can't understand their background supply chains.
# 
# First, we create a new database.

# In[20]:


aggregated_db = Database("EXIOBASE 2.2 aggregated")


# This is a normal database, not an `IOTable` database.

# In[21]:


type(aggregated_db)


# Now, we invert the EXIOBASE technosphere matrix.
# 
# This takes some minutes - around 4 on my laptop - so just be patient. It is helpful if there is plenty of free memory.

# In[22]:


inverse = np.linalg.pinv(dlca.technosphere_matrix.todense())


# With the inverse, we can calculated the aggregated inventories, and then write each aggregated process.

# In[23]:


inventory = dlca.biosphere_matrix * inverse
print(inventory.shape)


# Define the activity data fields we want to keep

# In[24]:


KEYS = (
    'exiobase_code',
    'group',
    'group_name',
    'location',
    'name',
    'synonym',
    'type',
    'unit'
)

data = {}


# Only take each non-zero biosphere flow, and create the aggregated processes.

# In[25]:


for ds in pyprind.prog_bar(Database("EXIOBASE 2.2")):
    col = dlca.activity_dict[ds.key]
    
    # Basic data
    data[("EXIOBASE 2.2 aggregated", ds['code'])] = {key: ds[key] for key in KEYS}
    # Exchanges
    data[("EXIOBASE 2.2 aggregated", ds['code'])]['exchanges'] = [{
        'type': 'biosphere',
        'amount': float(inventory[row, col]),
        'input': flow,
        'uncertainty type': 0
    } for flow, row in dlca.biosphere_dict.items() if inventory[row, col]]


# In[26]:


aggregated_db.write(data)


# We no longer need the dlca object, so we can forget about it to save some memory.

# In[27]:


dlca = None


# # Sample LCA calculations
# 
# We will look at two product systems selected in class. We found the dataset keys using code like:
# 
#     for x in Database("ecoinvent 2.2").search('fertili*'):
#         print(x, x.key)

# ## Cement production

# In[28]:


ex_cement = ('EXIOBASE 2.2 aggregated', 'Manufacture of cement, lime and plaster:CH')
ei_cement = ('ecoinvent 2.2', 'c2ff6ffd532415eda3eaf957b17b70a1')


# Check to make sure we have the correct activities

# In[31]:


get_activity(ex_cement)


# In[30]:


get_activity(ei_cement)


# In[32]:


lca = LCA({ex_cement: 1}, gwp)
lca.lci()
lca.lcia()
print("Exiobase:", lca.score / 1e6 / 10) # Assume 100 euros/ton

lca = LCA({ei_cement: 1}, gwp)
lca.lci()
lca.lcia()
print("Ecoinvent", lca.score)


# These numbers are remarkably similar. 
# 
# ## Nitrogenous fertilizer
# 
# Let's now look at nitrogen fertilizer:

# In[33]:


ei_n = ('ecoinvent 2.2', '920a20d9a87340557a31ee7e8a353d3c')
ex_n = ('EXIOBASE 2.2 aggregated', 'N-fertiliser:LU')


# Check to make sure we have the correct activities

# In[34]:


get_activity(ei_n)


# In[35]:


get_activity(ex_n)


# In[36]:


lca = LCA({ex_n: 1}, gwp)
lca.lci()
lca.lcia()
print("Exiobase:", lca.score / 1e6 * 0.8)  # Assume 800 euros/ton

lca = LCA({ei_n: 1}, gwp)
lca.lci()
lca.lcia()
print("Ecoinvent:", lca.score)


# This is quite interesting - more investigation would have to be done to understand why these values are so different.

# # Cleaning up
# 
# This project consumes a lot of hard drive space, about 2 gigabytes. We can get the exact size of this and all other projects (in gigabytes) with the following:

# In[37]:


projects.report()


# We can then delete the current project.
# 
# **This step is optional**, included as a convenience for those who do not want to work with Exiobase.

# In[38]:


projects.delete_project(delete_dir=True)


# The returned value is the name of the current project.

# In[39]:


projects.current

