#/usr/bin/bash
time mpirun.openmpi -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH --bynode -np 7 -hostfile hosts \
/shared/data/innocentive/dtra/bin/client cluster \
/shared/data/innocentive/nt.reduced.shuffled.fa /shared/data/innocentive/nt.reduced.shuffled20 \
/shared/data/innocentive/Testing/Testing1.fq /shared/data/innocentive/Testing/Testing1.lisa
