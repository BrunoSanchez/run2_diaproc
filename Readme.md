# Run2 DIA process 

## Steps to carry the analysis

* Create templates
    - Pick a tract-patch in Run2
    - Obtain the list of visits and filters for this patch
    - Configure the `coaddDriver.py`
    - Run the `coaddDriver.py`
    - Run the `multiBandDriver.py`
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
    - Beyond... 
