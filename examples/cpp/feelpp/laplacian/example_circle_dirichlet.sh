# the executable laplacian.e should exist

rootlorenzdir='../../../../'
executable='./'$rootlorenzdir'build/default/bin/laplacian.e'
configfile='circle/laplacian_dirichlet.cfg'
cmd=$executable' --config-file '$configfile
resultdir='/home/flecourtier/feelppdb/laplacian.e/np_1'

tab_h=(0.4 0.2 0.1 0.05 0.025)

echo '====== Execution with Feelpp ====='
declare -i i=0
for h in ${tab_h[*]}
do
    i+=1
    cmd_hsize=$cmd' --gmsh.hsize='$h
    echo '>>>> Execution with hsize='$h
    $cmd_hsize > circle/output/output_$i
    echo $cmd_hsize
    python3 write_values.py $i $h
done

echo '====== Convergence ====='
python3 convergence.py