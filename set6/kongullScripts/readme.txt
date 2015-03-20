How to compile and run on kongull cluster:

- Copy poisson.c, fst.f, CMakeLists.txt and all the jobscripts to kongull
- module load intelcomp
- module load openmpi/1.4.3-intel
- CC=icc FC=ifort cmake CMakeLists.txt
- make

Then we have compiled all we need and to run it we need to submit the jobs.
Run one of the following:
- qsub wrapper_n.sh
- qsub wrapper_p.sh
- qsub wrapper_t.sh

Output from each of the individal jobs come in a poisson_JC.xxx file where xxx specifies the job number. Gather output (ready to import to matlab) is aggregated to out.txt.
