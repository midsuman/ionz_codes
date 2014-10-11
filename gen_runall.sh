pbsfile="runall.pbs"

echo '#!/bin/bash' > $pbsfile
echo '#$ -N custom' >> $pbsfile
echo '#$ -cwd' >> $pbsfile
echo '#$ -pe openmpi 128' >> $pbsfile
echo '#$ -q mps.q' >> $pbsfile
echo '#$ -S /bin/bash' >> $pbsfile
# source modules environment:
echo "" >> $pbsfile
echo 'module add sge' >> $pbsfile
echo 'module add gcc/4.8.1' >> $pbsfile
echo 'module add intel-mpi/64/4.1.1/036' >> $pbsfile
echo 'module add gsl/gcc/1.15' >> $pbsfile

echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/c/cs/cs390/local/fftw-2.1.5/install/lib' >> $pbsfile

densdir="/research/prace/sph_smooth_cubepm_130315_6_1728_47Mpc_ext2/nc306/"
srcdir="/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/sources/"
cur_line=0
zlistfile="/mnt/lustre/scratch/cs390/47Mpc/snap_z3.txt"
z2listfile="/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
while read line
do
z3=$line
cur_line2=0

while read line2
do
if [ $cur_line = $cur_line2 ]; then
z2=$line2
fi
cur_line2=$((cur_line2+1))
done < $z2listfile

#echo "$z2   $z3"
echo 'echo z = $z3'  >> $pbsfile
echo 'mpirun -np $NSLOTS' "./ionz_main ${densdir}${z3}n_all.dat ${srcdir}${z2}.dat $z3" >> $pbsfile
#mpirun -np $NSLOTS ./ionz_main ${densdir}${z3}n_all.dat ${srcdir}${z2}.dat $z3
#mpirun -np $NSLOTS ../mpi_test
cur_line=$((cur_line+1))
done < $zlistfile
