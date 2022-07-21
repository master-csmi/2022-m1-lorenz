# the executable parareal.e should exist

# ./time_average.sh p
# p is the number of process

rootlorenzdir='../../../'
executable='./'$rootlorenzdir'build/default/bin/parareal.e'
# executable='mpirun -n 2 ./'$rootlorenzdir'build/default/bin/parareal.e'
filename='out.txt'

sum=0

for i in {1..10}
do
  $executable > $filename   
  n=1
  while read line; do
    # reading each line
    if [ $n -eq 2 ]
    then
        echo "Rep $i : $line"  
        sum=$((sum+line))
    fi
    n=$((n+1))
  done < $filename
done

echo "Sum : $sum" 
echo "Average (s) :"
bc -l <<< "$sum/10"

rm $filename