## command to run the multiBandDriver.py
## in order to make this possible we had to make soflinks to the
## actual data location in NERSC.
##
## The way to achieve this was, by making 

## `mkdir /global/cscratch1/sd/bos0109/run2.1i_softln/`
## `ln -s /global/cscratch1/sd/desc/DC2/data/Run2.1i/* /global/cscratch1/sd/bos0109/run2.1i_softln/.`
## `cp /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-dr4-v1/deepCoadd/skymap.pickle /global/cscratch1/sd/bos0109/templates_003/deepCoadd/.`

## edited the coadd repo yaml file so the parent directory is now in the
## directory we created before --> 
## _parents:
## - /global/cscratch1/sd/bos0109/run2.1i_softln/

## After that, we had to sele ct some options for the job to be able 
## to run, but happily it seems that it finally worked

nice -n 10 multiBandDriver.py /global/cscratch1/sd/bos0109/templates_003 --rerun multiband --id tract=4639 patch=0,0 filter=u^g^r^i^z^y  --job multiband  --cores 8 --time 5000  --batch-type=smp -C multibandDriver_config.py --config  measureCoaddSources.doPropagateFlags=False --clobber-config --batch-verbose  --batch-stats --batch-options='-C knl -q regular' --mpiexec='-bind-to socket'
