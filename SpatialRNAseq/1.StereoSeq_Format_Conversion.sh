############################################################
# Activate conda
############################################################

conda activate st


############################################################
# Read GEF files and convert into Seurat format (h5ad)
############################################################
#https://stereopy.readthedocs.io/en/latest/Tutorials/Format_Conversion.html

import stereo as st
import warnings
warnings.filterwarnings('ignore')


# read the GEF file
data_path = '/home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.gef'
data = st.io.read_gef(file_path=data_path, bin_size=20)
data.tl.cal_qc()
data.tl.raw_checkpoint()
# remember to set flavor as seurat
data = st.io.stereo_to_anndata(data,flavor='seurat',output='/home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.bin20.seurat_out.h5ad')


data_path = '/home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.gef'
data = st.io.read_gef(file_path=data_path, bin_size=1)
data.tl.cal_qc()
data.tl.raw_checkpoint()

# remember to set flavor as seurat
data = st.io.stereo_to_anndata(data,flavor='seurat',output='/home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.bin1.seurat_out.h5ad')



data_path = '/home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.gef'
data = st.io.read_gef(file_path=data_path, bin_size=200)
data.tl.cal_qc()
data.tl.raw_checkpoint()

# remember to set flavor as seurat
data = st.io.stereo_to_anndata(data,flavor='seurat',output='/home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.bin200.seurat_out.h5ad')



# read the GEF file
data_path = '/home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.gef'
data = st.io.read_gef(file_path=data_path, bin_size=100)
data.tl.cal_qc()
data.tl.raw_checkpoint()

# remember to set flavor as seurat
data = st.io.stereo_to_anndata(data,flavor='seurat',output='/home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.bin100.seurat_out.h5ad')



############################################################
# h5ad to RDS
############################################################

Rscript /home/woodydrylab/DiskArray/guest001/JW/bin/h5ad2rds.R \
	--infile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.bin20.seurat_out.h5ad \
	--outfile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/processed/JW/BGI_SAW_out/C04289G213.tissue.bin20.seurat_out.RDS


Rscript /home/woodydrylab/DiskArray/guest001/JW/bin/h5ad2rds.R \
	--infile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.bin1.seurat_out.h5ad \
	--outfile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/processed/JW/BGI_SAW_out/C04289G213.tissue.bin1.seurat_out.RDS


Rscript /home/woodydrylab/DiskArray/guest001/JW/bin/h5ad2rds.R \
	--infile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.bin200.seurat_out.h5ad \
	--outfile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/processed/JW/BGI_SAW_out/C04289G213.tissue.bin200.seurat_out.RDS


Rscript /home/woodydrylab/DiskArray/guest001/JW/bin/h5ad2rds.R \
	--infile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/visualization/C04289G213.bin200_1.0.h5ad \
	--outfile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/processed/JW/BGI_SAW_out/C04289G213.bin200_1.0.seurat_out.RDS


Rscript /home/woodydrylab/DiskArray/guest001/JW/bin/h5ad2rds.R \
	--infile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/outs/feature_expression/C04289G213.tissue.bin100.seurat_out.h5ad \
	--outfile /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/processed/JW/BGI_SAW_out/C04289G213.tissue.bin100.seurat_out.RDS




############################################################
# convert gene ID to Cluster
############################################################

#cd /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/processed/JW/BGI_SAW_out
cd /home/f06b22037/SSD2/JW/1136project_SingleCell/scRNA_stereo_visium
#ln -s /home/woodydrylab/FileShare/Stereoseq_raw_data/F24A040007172_PLAkisaT_0906/processed/JW/BGI_SAW_out/C04289G213.tissue.bin100.seurat_out.RDS ./
Rscript ./script/calculate_ortholog_UMI_counts_stereo.R all_group_long_convertedID_Four_Cla2_before20240211.csv C04289G213.tissue.bin100.seurat_out.RDS orthologUMI_StereoBin100_Cla.csv

