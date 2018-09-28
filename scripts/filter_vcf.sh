vcftools --vcf $1 --out out.1 --minDP 10 --minGQ 30 --recode --recode-INFO-all
vcftools --vcf out.1.recode.vcf --out out.2 --max-missing 0.9 --recode --recode-INFO-all
vcfallelicprimitives 
