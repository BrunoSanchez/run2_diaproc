nice -n 10 coaddDriver.py /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1 --output $SCRATCH/templates_002z --configfile $HOME/coadd_config_example.py  --id tract=4639 patch=0,0 filter=z --selectId visit=2340^2343^5890^6826^8007^8049^12448^12456^12484^32678^52543^159481^159494^167874^167875^169840^169848^174549^174602^177480^179970^179971^181868^181900^181938^181965^183892 --job templ_02z  --cores 32 --time 400  --batch-type=smp #  --batch-verbose  --batch-stats --batch-options='-C knl -q regular' --mpiexec='-bind-to socket' # --dry-run --clobber-output

nice -n 10 coaddDriver.py /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1 --output $SCRATCH/templates_002r --configfile $HOME/coadd_config_example.py  --id tract=4639 patch=0,0 filter=r --selectId visit=2340^2343^5890^6826^8007^8049^12448^12456^12484^32678^52543^159481^159494^167874^167875^169840^169848^174549^174602^177480^179970^179971^181868^181900^181938^181965^183892 --job templ_02r  --cores 32 --time 400  --batch-type=smp #  --batch-verbose  --batch-stats --batch-options='-C knl -q regular' --mpiexec='-bind-to socket' # --dry-run --clobber-input

nice -n 10 coaddDriver.py /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1 --output $SCRATCH/templates_002y --configfile $HOME/coadd_config_example.py  --id tract=4639 patch=0,0 filter=y --selectId visit=2340^2343^5890^6826^8007^8049^12448^12456^12484^32678^52543^159481^159494^167874^167875^169840^169848^174549^174602^177480^179970^179971^181868^181900^181938^181965^183892 --job templ_02y  --cores 32 --time 400  --batch-type=smp #  --batch-verbose  --batch-stats --batch-options='-C knl -q regular' --mpiexec='-bind-to socket' # --dry-run --clobber-output 

nice -n 10 coaddDriver.py /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1 --output $SCRATCH/templates_002g --configfile $HOME/coadd_config_example.py  --id tract=4639 patch=0,0 filter=g --selectId visit=2340^2343^5890^6826^8007^8049^12448^12456^12484^32678^52543^159481^159494^167874^167875^169840^169848^174549^174602^177480^179970^179971^181868^181900^181938^181965^183892 --job templ_02g  --cores 32 --time 400  --batch-type=smp #  --batch-verbose  --batch-stats --batch-options='-C knl -q regular' --mpiexec='-bind-to socket' # --dry-run --clobber-output 

nice -n 10 coaddDriver.py /global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1 --output $SCRATCH/templates_002i --configfile $HOME/coadd_config_example.py  --id tract=4639 patch=0,0 filter=i --selectId visit=2340^2343^5890^6826^8007^8049^12448^12456^12484^32678^52543^159481^159494^167874^167875^169840^169848^174549^174602^177480^179970^179971^181868^181900^181938^181965^183892 --job templ_02i  --cores 32 --time 400  --batch-type=smp #  --batch-verbose  --batch-stats --batch-options='-C knl -q regular' --mpiexec='-bind-to socket' # --dry-run --clobber-output 
