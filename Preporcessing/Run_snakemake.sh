
###############
# The command of running snakemake for running code from https://github.com/joweihsieh/Woodformation1136_SingleCell
###############


nohup snakemake --use-conda -c 48 TarCla2 > TarCla2_snakemake.log &
nohup snakemake --use-conda -c 48 Lchla2 > LchCla2_snakemake.log &
nohup snakemake --use-conda -c 48 EgrCla2 > EgrCla2_snakemake.log &
nohup snakemake --use-conda -c 48 PtrCla2 > PtrCla2_snakemake.log &

nohup snakemake --use-conda -c 48 Cla24 > Cla24_snakemake.log &
nohup snakemake --use-conda -c 48 Cla14 > Cla14_snakemake.log &
nohup snakemake --use-conda -c 48 Cla21 > Cla21_snakemake.log &


nohup snakemake --use-conda -c 48 PtrEgrTarLchCla > PtrEgrTarLchCla_snakemake.log &
