module load perftools-base
module load perftools

cd /work/ecse0803/ecse0803/es205/scaling/
cd cpl-library/
make clean
make
cd ../scaling
make clean
make md
rm -f md+samp
pat_build -o md+samp md

make cfd
rm -f cfd+samp
pat_build -o cfd+samp cfd

#qsub -q short qscript_perf
#aprun -n 96 ./cfd+samp : -n 96 ./md+samp

pat_report -o md_samp_1024.pat md+samp+*
pat_report -o cfd_samp_1024.pat cfd+samp+*


pat_report -O ca+src,load_balance -o md_furtherinfo_1024.pat md+samp+*
pat_report -O ca+src,load_balance -o cfd_furtherinfo_1024.pat cfd+samp+*
