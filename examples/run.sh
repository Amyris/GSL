alias gslc='c:/seq/GSL/gslc/bin/Debug/gslc'
alias gslc
export GSL_LIB='c:/seq/GSL/gslc/lib'
gslc --cm . simple_fusion --thumper simple_fusion --flat simple_fusion.txt --noprimers --ape . simple_fusion simple_fusion.gsl
gslc --cm . simple_megastitch --thumper simple_megastitch --flat simple_megastich.txt --ape . simple_megastitch simple_megastitch.gsl
gslc --cm . simple_promoter_gene_locus --thumper simple_promoter_gene_locus --flat simple_promoter_gene_locus.txt --ape . simple_promoter_gene_locus simple_promoter_gene_locus.gsl
gslc --cm . terpene_design --thumper terpene_design --flat terpene_design.txt --ape . terpene_design terpene_design.gsl
gslc --cm . allele_swap --thumper allele_swap --flat allele_swap.txt --ape . allele_swap allele_swap.gsl


#gslc --flat simple_fusion.txt --noprimers --ape . simple_fusion simple_fusion.gsl
#gslc --flat simple_megastich.txt --ape . simple_megastitch simple_megastitch.gsl
#gslc --flat simple_promoter_gene_locus.txt --ape . simple_promoter_gene_locus simple_promoter_gene_locus.gsl
#gslc --flat terpene_design.txt --ape . terpene_design terpene_design.gsl
#gslc --flat allele_swap.txt --ape . allele_swap allele_swap.gsl
