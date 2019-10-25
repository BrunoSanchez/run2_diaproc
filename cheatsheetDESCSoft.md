A cheat sheet for lsst-desc software APIs
-----------------------------------------

## GCR:  General catalog reader tool

This grants access to datasets in catalog format.  
It performs several tasks such as querying and filtering.

### The way of...

* importing it: 
`import GCRCatalogs`

* List every available catalog
`all_catalogs = list(GCRCatalogs.available_catalogs.keys())`

* load cats with it
`cat = GCRCatalogs.load_catalog('dc2_truth_run1.2_static')`

* list their columns
`cat.list_all_quantities(include_native=True)`

* Find info on a column
`info_dict = cat.get_quantity_info(col_name)`

* Actually load data
`cat_data = cat.get_quantities(['qty1', 'qty2'], 
                               native_filters=['qty3>0'],
                               filters=[(lambda x:x=='bla', 'qty4')])`
* loading as an iterator
`return_iterator=True` flag


### Native qtys and filters

It seems that the already existent quantities in the catalogs are much easier to use 
since they are faster for querying and filtering. Any other non-native qty seems to be
calculated on the fly.

## Butler: the tool for data repos access

This allows you to retrieve data from the stack, and is the preferred way of interaction
with the datasets according to **dm**.

### The way of...

* importing the module
`from lsst.daf.persistence import Butler`

* instance a new butler
`butler = Butler('/path/to/data/repo')`

* list the data available in the butler
    First need the _mapper_ class that the butler is using:  
     `mapper = butler._getDefaultMapper()(repo)`
    then you can actually look into the dataset types available
     `all_diadataset_types = diamapper.getDatasetTypes()`
    
* getting images


* displaying images

### Available things to work with

* The skymap: `skymap = butler.get('deepCoadd_skyMap')`
This object contains several ways to access the corresponding tract and patches that build up the coadd.  
For instance one can use the `skymap.findTractPatchList` method, in order to obtain the overlapping t+p given a set of coordinates. The set of coordinates should be a list of `SpherePoint` objects, with the ICRS coordinates.