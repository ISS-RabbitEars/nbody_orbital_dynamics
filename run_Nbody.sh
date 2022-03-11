
var=1
while [ $var -eq 1 ]
do
python3.9 gravity_Nbody_OD_ode.py && var=0 || var=1
done
