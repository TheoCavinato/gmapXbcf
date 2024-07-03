
VCF=/home/theo/Desktop/Data/GSA_filtered_1000GP_va/GSA_chr22.bcf
MAP=~/commands/GLIMPSE/maps/genetic_maps.b37/chr22.b37.gmap.gz
OUT=new_map.plink_format
python3 gmapXbcf.py --vcf $VCF --map $MAP --out $OUT --format plink