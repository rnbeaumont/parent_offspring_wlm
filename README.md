# Estimate independent Parent-Offspring Genetic Associations with arbitrary sample overlap

### Usage:
perl apply_wlm.pl \
  --child_file <path to file containing child association estimates> \
  --mother_file <path to file containing mother association estimates> \
  --father_file <path to file containing father association estimates> (only required if allowing non-zero paternal effects) \
  --cm <child/mother overlap> \
  --cp <child/father overlap> (only required if allowing non-zero paternal effects) \
  --mp <mother/father overlap> (only required if allowing non-zero paternal effects) \
  --out_file <output filename>
  
  where the child, mother and father files should contain columns (in order):
    marker name
    effect allele
    other allele
    beta
    se
    p
  and the overlap estimates provided to --cm, --cp and --mp should be the proportion of individuals overlapping between each pair of GWAS multiplied by the phenotypic correlation. This term can be estimated as the intercept of an LD-score regression analysis.
