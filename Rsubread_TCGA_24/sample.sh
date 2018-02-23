#!/usr/bash
awk '{if(NR==1)print $0}' GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt>header
for ((i=1;i<=9264;i++))
do
awk '{print $"'$i'"}' header
done
