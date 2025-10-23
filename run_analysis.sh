#!/bin/bash
#

analysis=poisson

for i in divide50p # full sc25 sc10 cc25 cc10 divide10 divide10p divide20 divide20p divide50 divide50p
do
  if [ $i == full ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson.sh
      sed -i "s/memory/600GB/g" run_poisson.sh
      sbatch run_poisson.sh
      sed -i "s/$i/option/g" run_poisson.sh
      sed -i "s/600GB/memory/g" run_poisson.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif.sh
      sed -i "s/memory/650GB/g" run_cif.sh
      sbatch run_cif.sh
      sed -i "s/$i/option/g" run_cif.sh
      sed -i "s/650GB/memory/g" run_cif.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/75GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/75GB/memory/g" run_cox.sh

    fi

  elif [ $i == sc25 ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson.sh
      sed -i "s/memory/600GB/g" run_poisson.sh
      sbatch run_poisson.sh
      sed -i "s/$i/option/g" run_poisson.sh
      sed -i "s/600GB/memory/g" run_poisson.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif.sh
      sed -i "s/memory/300GB/g" run_cif.sh
      sbatch run_cif.sh
      sed -i "s/$i/option/g" run_cif.sh
      sed -i "s/300GB/memory/g" run_cif.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/30GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/30GB/memory/g" run_cox.sh

    fi
    
 elif [ $i == sc10 ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson.sh
      sed -i "s/memory/300GB/g" run_poisson.sh
      sbatch run_poisson.sh
      sed -i "s/$i/option/g" run_poisson.sh
      sed -i "s/300GB/memory/g" run_poisson.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif.sh
      sed -i "s/memory/60GB/g" run_cif.sh
      sbatch run_cif.sh
      sed -i "s/$i/option/g" run_cif.sh
      sed -i "s/60GB/memory/g" run_cif.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/30GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/30GB/memory/g" run_cox.sh

    fi

  elif [ $i == cc25 ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson.sh
      sed -i "s/memory/600GB/g" run_poisson.sh
      sbatch run_poisson.sh
      sed -i "s/$i/option/g" run_poisson.sh
      sed -i "s/600GB/memory/g" run_poisson.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif.sh
      sed -i "s/memory/400GB/g" run_cif.sh
      sbatch run_cif.sh
      sed -i "s/$i/option/g" run_cif.sh
      sed -i "s/400GB/memory/g" run_cif.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/30GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/30GB/memory/g" run_cox.sh

    fi

  elif [ $i == cc10 ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson.sh
      sed -i "s/memory/400GB/g" run_poisson.sh
      sbatch run_poisson.sh
      sed -i "s/$i/option/g" run_poisson.sh
      sed -i "s/400GB/memory/g" run_poisson.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif.sh
      sed -i "s/memory/300GB/g" run_cif.sh
      sbatch run_cif.sh
      sed -i "s/$i/option/g" run_cif.sh
      sed -i "s/300GB/memory/g" run_cif.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/30GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/30GB/memory/g" run_cox.sh

    fi

  elif [ $i == divide10 ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson-divide.sh
      sed -i "s/memory/600GB/g" run_poisson-divide.sh
      sbatch run_poisson-divide.sh
      sed -i "s/$i/option/g" run_poisson-divide.sh
      sed -i "s/600GB/memory/g" run_poisson-divide.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif-divide.sh
      sed -i "s/memory/50GB/g" run_cif-divide.sh
      sbatch run_cif-divide.sh
      sed -i "s/$i/option/g" run_cif-divide.sh
      sed -i "s/50GB/memory/g" run_cif-divide.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/30GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/30GB/memory/g" run_cox.sh

    fi

    elif [ $i == divide20 ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson-divide.sh
      sed -i "s/memory/600GB/g" run_poisson-divide.sh
      sbatch run_poisson-divide.sh
      sed -i "s/$i/option/g" run_poisson-divide.sh
      sed -i "s/600GB/memory/g" run_poisson-divide.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif-divide.sh
      sed -i "s/memory/50GB/g" run_cif-divide.sh
      sbatch run_cif-divide.sh
      sed -i "s/$i/option/g" run_cif-divide.sh
      sed -i "s/50GB/memory/g" run_cif-divide.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/30GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/30GB/memory/g" run_cox.sh

    fi

    elif [ $i == divide50 ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson-divide.sh
      sed -i "s/memory/600GB/g" run_poisson-divide.sh
      sbatch run_poisson-divide.sh
      sed -i "s/$i/option/g" run_poisson-divide.sh
      sed -i "s/600GB/memory/g" run_poisson-divide.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif-divide.sh
      sed -i "s/memory/50GB/g" run_cif-divide.sh
      sbatch run_cif-divide.sh
      sed -i "s/$i/option/g" run_cif-divide.sh
      sed -i "s/50GB/memory/g" run_cif-divide.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/30GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/30GB/memory/g" run_cox.sh

    fi

   elif [ $i == divide10p ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson-divide.sh
      sed -i "s/memory/200GB/g" run_poisson-divide.sh
      sbatch run_poisson-divide.sh
      sed -i "s/$i/option/g" run_poisson-divide.sh
      sed -i "s/200GB/memory/g" run_poisson-divide.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif-divide.sh
      sed -i "s/memory/200GB/g" run_cif-divide.sh
      sbatch run_cif-divide.sh
      sed -i "s/$i/option/g" run_cif-divide.sh
      sed -i "s/200GB/memory/g" run_cif-divide.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/200GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/200GB/memory/g" run_cox.sh

    fi

elif [ $i == divide20p ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson-divide.sh
      sed -i "s/memory/100GB/g" run_poisson-divide.sh
      sbatch run_poisson-divide.sh
      sed -i "s/$i/option/g" run_poisson-divide.sh
      sed -i "s/100GB/memory/g" run_poisson-divide.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif-divide.sh
      sed -i "s/memory/300GB/g" run_cif-divide.sh
      sbatch run_cif-divide.sh
      sed -i "s/$i/option/g" run_cif-divide.sh
      sed -i "s/300GB/memory/g" run_cif-divide.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/200GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/200GB/memory/g" run_cox.sh

    fi

 elif [ $i == divide50p ]
  then

    if [ $analysis == poisson ]
    then
      sed -i "s/option/$i/g" run_poisson-divide.sh
      sed -i "s/memory/300GB/g" run_poisson-divide.sh
      sbatch run_poisson-divide.sh
      sed -i "s/$i/option/g" run_poisson-divide.sh
      sed -i "s/300GB/memory/g" run_poisson-divide.sh

    elif [ $analysis == cif ]
    then
      sed -i "s/option/$i/g" run_cif-divide.sh
      sed -i "s/memory/600GB/g" run_cif-divide.sh
      sbatch run_cif-divide.sh
      sed -i "s/$i/option/g" run_cif-divide.sh
      sed -i "s/600GB/memory/g" run_cif-divide.sh

    elif [ $analysis == cox ]
    then
      sed -i "s/option/$i/g" run_cox.sh
      sed -i "s/memory/200GB/g" run_cox.sh
      sbatch run_cox.sh
      sed -i "s/$i/option/g" run_cox.sh
      sed -i "s/200GB/memory/g" run_cox.sh

    fi

  fi

done
