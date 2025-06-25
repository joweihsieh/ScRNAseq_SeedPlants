repair.sh -Xmx20g in1=SRR28089563_S1_L001_R1_001.fastq in2=SRR28089563_S1_L001_R2_001.fastq out1=R1_fixed.fastq out2=R2_fixed.fastq outs=singles.fastq repair

mv R1_fixed.fastq SRR28089563_S1_L001_R1_001.fastq
mv R2_fixed.fastq SRR28089563_S1_L001_R2_001.fastq



spaceranger count  \
      --id="spaceranger_count_TenX_stemQuanzi2024_v3" \
      --transcriptome=/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/spaceranger_ref_TenX_Ptr \
      --fastqs=/home/woodydrylab/FileShare/Spatialtranscriptomes/fastq3   \
      --sample=SRR28089563 \
      --localcores=10 \
      --create-bam=true \
      --localmem=128 \
      --localvmem=128 \
      --slide=V12A25-388 \
      --area=C1 \
      --image=/home/woodydrylab/FileShare/10Xvisium_Quanzi/image/WT-2-9.jpg

