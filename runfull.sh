# module load pe_archive
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

cd dia_pipe
scons
cd ..

setup -jr dia_pipe/

git clone https://github.com/lsst/obs_lsst.git
cd obs_lsst/
git checkout w.2019.19_diff
scons
cd ..
setup -jr obs_lsst/

export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1

#imageDifferenceDriver.py /global/cscratch1/sd/desc/DC2/data/Run1.2p/w_2018_39/rerun/coadd-v4/ --output test_imdiff001  --id     visit=431306  -C dia_pipe/config/imageDifferenceDriver.py --batch-type=slurm --mpiexec='-bind-to socket'   --cores 10 --job test_001 --time 500 --batch-options='-C knl -q regular'

#imageDifferenceDriver.py /global/cscratch1/sd/desc/DC2/data/Run1.2p/w_2018_39/rerun/coadd-v4/ --output test_imdiff002  --id     visit=431306  -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=slurm --mpiexec='-bind-to socket'   --cores 10 --job test_002 --time 500 --batch-options='-C knl -q regular'


#sqlite3 /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1/tracts_mapping.sqlite3
#select distinct(visit) from overlaps where tract=4849 and patch='(5, 5)' order by visit;
#  output
# 190265
# 227032
# 399393
# 458495
# 512050
# 680296
# 687270


################################################################################
############# Full set of visits
################################################################################
# with ALupton
time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_full  --id     visit=190265^227032^399393^458495^512050^680296^687270  -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'

time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_visit190265  --id     visit=190265  -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'

time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_visit227032 --id visit=227032 -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'

time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_visit399393 --id visit=399393 -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_visit458495 --id visit=458495  -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_visit512050 --id visit=512050 -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_visit680296 --id visit=680296 -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket' --cores 32 --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_visit687270  --id     visit=687270  -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


# with ZOGY
time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_fullZOGY  --id     visit=190265^227032^399393^458495^512050^680296^687270  -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'

time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_ZOGYvisit190265  --id     visit=190265  -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'

time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_ZOGYvisit227032 --id visit=227032 -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_ZOGYvisit399393 --id visit=399393 -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_ZOGYvisit458495 --id visit=458495  -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_ZOGYvisit512050 --id visit=512050 -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_ZOGYvisit680296  --id     visit=680296  -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


time imageDifferenceDriver.py  /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/coadd-v1/ --output /global/cscratch1/sd/bos0109/test_imdiff_ZOGYvisit687270  --id     visit=687270  -C dia_pipe/config/imageDifferenceDriverZOGY.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'


