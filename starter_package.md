
## Starter package for DIA subtractions using the LSST stack.

Hi! This file is to share some of the details regarding DIA analysis in NERSC/Cori using the lsst pipeline stack provided by DESC. 
In principle, you don't need to know much of the configuration 
----------

### Environment setup 
#### Pipeline DIA stack
First of all you will need to set up the environment. 
From now on I am assuming you have loaded the needed python environments that are listed in the confluence page dedicated to starting on NERSC and jupyter. If you haven't done this, please come back after this is completed.

The shell command to run is the following, and it needs not to be run from 
any specific directory:

`source /global/cfs/cdirs/lsst/groups/SN/dia/code/v19/setup_nersc.sh`

After this you should be able to see the prefix of your prompt with the tag of the 
pipeline build you are using.

It is rather convienent to set this as an alias in your `.bashrc` file in your home directory at Cori.

#### DIA Shared area at NERSC filesystem
The computing infrastructure DESC team has prepared a specific area in the filesystem of Cori to place data and code at the following locations:
`/global/cfs/cdirs/lsst/groups/SN/dia`

It is also important to know, that the previous script sets up a couple of environment variables that hold this handy information for you:
```
DESC_DIA_DIR = /global/cfs/cdirs/lsst/groups/SN/dia/code/v19
DESC_DIA_INSTALL = /global/cfs/cdirs/lsst/groups/SN/dia/code
DESC_DIA_VER = v19
```
It is also recommended that you create your own directory for storing output of your analysis under the following location:
`$DESC_DIA_DIR/../../data/`

At this locations the permissions are open for any member of the collaboration. This makes sharing results a lot easier and collaborating inside the DIA topical team.

### The templates

I have created a bunch of templates using Run2.1i data, for the rectangular region of `56 < RA < 58; -32 < Dec < -31` and it can be used anytime for testing subtraction code.

The templates are located in this path
`/global/cscratch1/sd/bos0109/templates_rect`

As we can see the directory structure looks like this

```
templates_rect
    ... config  
    ... deepCoadd  
    ... deepCoadd-results  
    ... repositoryCfg.yaml  
    ... rerun  
    ... schema
```
The template images are located in the `deepCoadd` subdirectory.   
The data is splitted into filters and after that the so called _tracts_ and _patches_, 
which are the sky tiling that LSST is going to use.

To check the image for instance for tract=4431 and patch=0,3 filter=g you need to look into this file:

`ds9 /global/cscratch1/sd/bos0109/templates_rect/deepCoadd/g/4431/0,3.fits`

### The subtraction command

In order to make subtractions the command is the following 

```
imageDifferenceDriver.py  [template_dir] --rerun [output_dirname] --id visit=[visitnumber] detector=[0...179] -C ./config/imageDifferenceDriver_config.py --batch-type=smp --cores 8 --job [job_name] --time 100  
```

The argument `--batch-type=smp` means you are running at the login node in cori. If, instead, `--batch-type=slurm` you are going to run in the cori batch queue system, and you need to add this flags `--batch-verbose  --batch-stats --mpiexec='-bind-to socket' --batch-options='-C haswell -q shared'`. The last flag, actually, selects the cpus to use, cori has `knl` nodes as well as `haswell` nodes. Haswell is more expensive, and both of them have `shared` or `regular` queues (the after `-q ` flag).

In summary:
* for haswell, shared queue use `--batch-options='-C haswell -q shared`
* for knl, regular queue `--batch-options='-C knl -q regular`

An example command, that takes 12 minutes to run is: 

```
time nice -n 5 imageDifferenceDriver.py  /global/cscratch1/sd/bos0109/templates_rect --rerun diff_test --id visit=193822 -C ./config/imageDifferenceDriver_config.py --batch-type=smp --cores 8 --job imdifftest_v193822_fg --time 100   --batch-verbose  --batch-stats --mpiexec='-bind-to socket'
```


### The output files

In every subtraction, under the results output directory we shall find something like this

```
diff_test
    ... config
    ... deepDiff
    ... repositoryCfg.yaml 
    ... schema
```
and the images, and individual detection tables are in
```
deepDiff
    ... v00193822-fg
        ... R30
            ... diaSrc_00193822-g-R30-S10-det120.fits 
            ... diffexp_00193822-g-R30-S10-det120.fits
            ... diaSrc_00193822-g-R30-S11-det121.fits 
            ... diffexp_00193822-g-R30-S11-det121.fits
            ... diaSrc_00193822-g-R30-S12-det122.fits  
            ... diffexp_00193822-g-R30-S12-det122.fits
            ... diaSrc_00193822-g-R30-S20-det123.fits  
            ... diffexp_00193822-g-R30-S20-det123.fits
            ... diaSrc_00193822-g-R30-S21-det124.fits 
            ... diffexp_00193822-g-R30-S21-det124.fits
            ... diaSrc_00193822-g-R30-S22-det125.fits 
            ... diffexp_00193822-g-R30-S22-det125.fits
        ... R31
        ... R32
        ... R41
        ... R42
```

