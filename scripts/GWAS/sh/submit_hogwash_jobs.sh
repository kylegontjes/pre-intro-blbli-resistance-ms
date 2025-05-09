#!/bin/sh

cd /nfs/turbo/umms-esnitkin/Project_Penn_KPC/Analysis/Combination_resistance/pre-intro-blbli-resistance-ms/scripts/GWAS/sh

for script in `ls *hogwash*.sbat | grep -v model_hogwash  `
  do
  echo $script
  sbatch $script
done

