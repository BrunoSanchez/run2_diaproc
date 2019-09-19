## dia commands

time imageDifferenceDriver.py  /global/cscratch1/sd/bos0109/templates_003/rerun/ --output /global/cscratch1/sd/bos0109/test_imdiff_run2  --id visit=  -C dia_pipe/config/imageDifferenceDriver.py --batch-type=smp --mpiexec='-bind-to socket'   --cores 32  --job test_full --time 5000 --batch-options='-C knl -q regular'
