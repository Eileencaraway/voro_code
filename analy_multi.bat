./compile.x
for((i=1;i<10;i+=1))
do 
./a.out 'dumpfile.00000'$i'0000.txt'
 
mv local_density.dat 'local_density-'$i'0000.dat' 
done 
