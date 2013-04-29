#/usr/bin/bash
time ~/install/bin/mpirun -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH --bynode -np 7 -hostfile hosts /home/fpavetic/dtra/bin/client cluster /home/fpavetic/db/nt10000 /home/fpavetic/dtra/index10000 /home/fpavetic/db/nt10000_500 /home/fpavetic/dtra/result
