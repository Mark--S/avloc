#$$ -S /bin/bash
#$$ -o /epp/scratch/neutrino/ms711/logs
#$$ -e /epp/scratch/neutrino/ms711/logs
#$$ -q mps.q
source /mnt/lustre/epp_scratch/neutrino/rat/snoing/install/env_rat-5.0.0.sh
rat -l ${Log} -o ${Output} ${Macro}
