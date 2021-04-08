# parent_offspring_wlm

## script to apply WLM to parent/offspring data

### Usage:
perl apply_wlm.pl \
  --child_file <path to file containing child association estimates> \
  --mother_file <path to file containing mother association estimates> \
  --father_file <path to file containing father association estimates> \
  --cm <child/mother overlap> \
  --cp <child/father overlap> \
  --mp <mother/father overlap> \
  --out_file <output filename>
  
  where the child, mother and father files should contain columns:
    marker name
    effect allele
    other allele
    beta
    se
    p
  and the overlap estimates provided to --cm, --cp and --mp should be the proportion of individuals overlapping between each pair of GWAS multiplied by the phenotypic correlation. This term can be estimated as the intercept of an LD-score regression analysis
