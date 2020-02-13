
## Starter package for DIA subtractions using the LSST stack.

----------

### Environment setup 
First of all you will need to set up the environment. 
The set of shell commands to do this are the following, and they need to run from 
your home directory at NERSC CORI computers:

```
module unload python
module unload PrgEnv-intel/6.0.5
module load PrgEnv-gnu/6.0.5
module swap gcc gcc/8.3.0
module rm craype-network-aries
module rm cray-libsci/19.02.1
module unload craype
export CC=gcc
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2019_19/loadLSST.bash
setup lsst_distrib
cd $HOME/dia_pipe
scons
cd ..
setup -jr $HOME/dia_pipe/
git clone https://github.com/lsst/obs_lsst.git
cd $HOME/obs_lsst/
git checkout w.2019.19_diff
scons
cd ..
setup -jr $HOME/obs_lsst/
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
```

After this you should be able to see the prefix of your prompt with the tag of the 
pipeline build you are using.

Once you do this once, the lines with the `git clone` and `scons` they can be ignored, 
saving setup time for your following terminal sessions.
The setup script then should read:
```
module unload python
module unload PrgEnv-intel/6.0.5
module load PrgEnv-gnu/6.0.5
module swap gcc gcc/8.3.0
module rm craype-network-aries
module rm cray-libsci/19.02.1
module unload craype
export CC=gcc
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2019_19/loadLSST.bash
setup lsst_distrib
#cd $HOME/dia_pipe
#scons
#cd ..
setup -jr $HOME/dia_pipe/
#git clone https://github.com/lsst/obs_lsst.git
#cd $HOME/obs_lsst/
#git checkout w.2019.19_diff
#scons
#cd ..
setup -jr $HOME/obs_lsst/
export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1
```

### The templates

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

