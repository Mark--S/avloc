#$$ -S /bin/bash
#$$ -o /epp/scratch/neutrino/sjp39/logs
#$$ -e /epp/scratch/neutrino/sjp39/logs
#$$ -q mps.q
cd /epp/scratch/neutrino/sjp39
source env_rat-dev.sh
rat -l ${Log} -o ${Output} ${Macro}
