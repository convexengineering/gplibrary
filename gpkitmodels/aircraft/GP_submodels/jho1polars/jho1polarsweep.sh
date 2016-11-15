AIRFOIL=jho1
Re="200 250 300 350 400 450 500 550 600 650 700"
for r in $Re
do
    ./genpolar.sh $AIRFOIL $r
done
