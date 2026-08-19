#!/bin/bash;

echo -e "\n";
echo "######################################";
echo "#Variant Annotation and Filteration Pipeline#";
echo "######################################";

input_vcf = "${1}" 2> warnings.txt

#Root directories
pipeline_root="/mnt/g/Projects/Divya/Final/";
snpeff_snpsift_root="/mnt/g/Projects/Divya/snpEff/";

echo -e "\n";
echo "###################################";
echo "#Filtering variants on the quality parameters#";
echo "###################################";

java -Xmx8G -jar "${snpeff_snpsift_root}"SnpSift.jar filter "(FILTER = 'PASS')" "${1}" > "${1}"_Quality_filtered.vcf

echo "Filtering step done!"

echo -e "\n";
echo "###################################";
echo "#Adding dbSNP Ids to variants#";
echo "###################################";

java -Xmx8G -jar "${snpeff_snpsift_root}"SnpSift.jar annotate -id dbsnp.vcf.gz "${1}"_Quality_filtered.vcf > "${1}"_QF_dbsnp.vcf 2> warnings.txt

echo "dbSNP ID addition step done!"

echo -e "\n";
echo "###################################";
echo "#Annotating with SnpEff on GRCh38#";
echo "###################################";

java -Xmx8G -jar "${snpeff_snpsift_root}"snpEff.jar GRCh38.p13.RefSeq "${1}"_QF_dbsnp.vcf > "${1}"_QF_dbsnp_ann.vcf


echo "Functional annotation step done!"

echo -e "\n";
echo "###################################";
echo "#Annotating with ClinVar#";
echo "###################################";

java -Xmx8G -jar "${snpeff_snpsift_root}"SnpSift.jar annotate clinvar_chr_corrected_latest.vcf.gz "${1}"_QF_dbsnp_ann.vcf > "${1}"_QF_dbsnp_ann_clinvar.vcf

echo "ClinVar annotation step done!"

echo -e "\n";
echo "###################################";
echo "#Annotating variants with GWAS Catalog#";
echo "###################################";

java -Xmx8G -jar "${snpeff_snpsift_root}"SnpSift.jar gwasCat -db gwas_catalog_v1.0-associations_e100_r2021-06-08.tsv "${1}"_QF_dbsnp_ann_clinvar.vcf > "${1}"_QF_dbsnp_ann_clinvar_gwas.vcf

echo "GWAS Annotation done!"

echo -e "\n";
echo "###################################";
echo "#Post Annotation filtering#";
echo "###################################";

#Variants only present in dbSNP
java -Xmx8G -jar "${snpeff_snpsift_root}"SnpSift.jar filter -f "${1}"_QF_dbsnp_ann_clinvar_gwas.vcf "exists ID" > "${1}"_ann_variants_in_dbsnp.vcf

#Variant with MAF > 0.01

#java -jar "${snpeff_snpsift_root}"SnpSift.jar filter "(AF_TGP >= 0.01)" "${1}"_QF_ann_dbsnp_clinvar_in_dbsnp.vcf > "${1}"_QF_ann_dbsnp_clinvar_in_dbsnp_AF.vcf

echo "Post-annotation filtering step done!"

echo -e "\n";
echo "###################################";
echo "#Creating a table with required fields#";
echo "###################################";

cat "${1}"_ann_variants_in_dbsnp.vcf | "${snpeff_snpsift_root}"scripts/vcfEffOnePerLine.pl | java -jar "${snpeff_snpsift_root}"SnpSift.jar extractFields - CHROM ID POS REF ALT "GEN[*].GT" "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].HGVS_C" CLNSIG CLNVC CLNDN CLNREVSTAT AF_TGP > "${1}"_all_variants_table.txt

cat "${1}"_ann_variants_in_dbsnp.vcf | "${snpeff_snpsift_root}"scripts/vcfEffOnePerLine.pl | java -jar "${snpeff_snpsift_root}"SnpSift.jar extractFields - CHROM ID POS REF ALT "ANN[*].GENE" GWASCAT_TRAIT GWASCAT_P_VALUE GWASCAT_OR_BETA > "${1}"_all_gwas_table.txt

echo "Variant table creation done!"

echo -e "\n";
echo "###################################";
echo "#Manupulating the table of variants to get variants of interest#";
echo "###################################";

#remove duplicates
sort "${1}"_all_variants_table.txt | uniq > "${1}"_uniq_variants.txt
sort "${1}"_all_gwas_table.txt | uniq > "${1}"_all_gwas_variants.txt

#Replace genotype values from numbers (e.g. 1/1) to letters (e.g. AA)

bcftools query -f '%ID\t[%TGT]\n' "${1}"_ann_variants_in_dbsnp.vcf > gt.txt

awk -vOFS="\t" 'NR==FNR{a[$1]=$2; next}{$6=a[$2]; print}' gt.txt "${1}"_uniq_variants.txt > gt_replaced.txt

awk 'NR==1{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\tGT\tEFFECT\tGENE\tIMAPCT\tHGVS_C\t"$10"\t"$11"\t"$12"\t"$13"\t"$14;next}{print}' gt_replaced.txt > "${1}"_actionable_variants.txt

#get the file header
awk 'NR==1{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\tGT\tEFFECT\tGENE\tIMAPCT\tHGVS_C\t"$10"\t"$11"\t"$12"\t"$13"\t"$14;next}' gt_replaced.txt > header

#Filter on the genelist of selected conditions

for i in `awk '{print $1}' Genelist.txt`; do grep "$i" -w "${1}"_actionable_variants.txt; done > "${1}"_variants_genelist.txt

cat header "${1}"_variants_genelist.txt > newfile
mv newfile "${1}"_variants_genelist.txt

#all pathogenic
cat "${1}"_actionable_variants.txt | grep -E "Pathogenic|Likely_pathogenic" > "${1}"_variants_pathogenic.txt
cat header "${1}"_variants_pathogenic.txt > newfile
mv newfile "${1}"_variants_pathogenic.txt

#all risk factor
cat "${1}"_actionable_variants.txt | grep "risk_factor" > "${1}"_variants_riskfactor.txt
cat header "${1}"_variants_riskfactor.txt > newfile
mv newfile "${1}"_variants_riskfactor.txt

#all protective
cat "${1}"_actionable_variants.txt | grep "protective" > "${1}"_variants_protective.txt
cat header "${1}"_variants_protective.txt > newfile
mv newfile "${1}"_variants_protective.txt

#all drug response
cat "${1}"_actionable_variants.txt | grep "drug_response" > "${1}"_variants_drug_response.txt
cat header "${1}"_variants_drug_response.txt > newfile
mv newfile "${1}"_variants_drug_response.txt

rm warnings.txt
#rm gt.txt
#rm gt_replaced.txt
#rm header
rm "${1}"_all_variants_table.txt
rm "${1}"_uniq_variants.txt
rm "${1}"_all_gwas_table.txt

echo "Variant Annotation Pipeline finished!"
