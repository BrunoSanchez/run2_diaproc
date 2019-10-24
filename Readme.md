# Run2 DIA process 
This repo contains several aids for script creation, that make use of the _drivers_ available in `dia_pipe` repo.  
The intention is to have a guide so anyone can build template images from calexps, and after that perform DIA analysis.  
The scripts provided work as a _cli_, and they would ingest a tract and a patch numbers, and produce _shell_ files with the commands so the user can run this at NERSC.

Several options are available, but the user needs to specify the following information:
* `database`: the `tracts_mapping.sqlite3` SQLite database for tract/patch selection
* `calexp_repo`: the `calexp` butler repository, so we can use the images there to produce template coadds.
* `tract, patch`
* Depending on the step of the analysis you will need to set the `templ` or template repo, and the `dia_repo`


## Steps to carry the analysis

* Create templates
    - Pick a tract-patch in Run2
    - Obtain the list of visits and filters for this patch
    - Configure the `coaddDriver.py`
    - Run the `coaddDriver.py`
    - Run the `multiBandDriver.py` on the coadd repo
* Run the DIA pipeline
    - Configure the chosen DIA algorithm 
    - Run `imageDifferenceDriver.py`
    - Obtain the `DIASources` and `DIAObjects` by means of `associationDriver.py`
    - Filter the detections to create sensible catalogs 
    - Optionally create stamps, lightcurves and ...
* Study results
    - Obtain a matching in space-time between the `DIAObjects` and truth catalogs
    - Obtain metric figures quantifying the completitude-contamination 
    - Measure selection functions
    - Etc.
