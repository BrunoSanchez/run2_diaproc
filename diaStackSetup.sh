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
