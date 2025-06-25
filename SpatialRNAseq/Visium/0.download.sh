#mkdir -p /home/woodydrylab/FileShare/Spatialtranscriptomes
#cd /home/woodydrylab/FileShare/Spatialtranscriptomes

wget -q --show-progress https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR28089563/SRR28089563 -O SRR28089563.sra
vdb-validate SRR28089563.sra
fastq-dump --split-files --origfmt --gzip SRR28089563.sra

mv SRR28089563_1.fastq SRR28089563_S1_L001_R1_001.fastq
mv SRR28089563_2.fastq SRR28089563_S1_L001_R2_001.fastq



######
######
###### installation
mkdir spaceranger
cd spaceranger

curl -o spaceranger-3.1.1.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.1.tar.gz?Expires=1726585324&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=btiTmqqrXZ6sSB~luRk8TMfaE5y3hsslgPVkHnRrkeVuk-E8b58tsG87Wr3m642Vd2z-L74WQz6M4uQ7yE85lsDA6iLiMCcYggO3IiDF9fzR3tzxC7x3BjCEQY6nsC57tkQcpqvR7XSlNXZ4BRL626eKNr5pfp79pR4dB25sdrHPj628PlkUPIcMPIBct5YAgzz-WVWO~RamkvMu57F5HdMphetw4ZPYGXK3DSKZhq453YmaDX8RISrND6NZSJdhdGBYpGCtOutnlHKTk4tF2GurX7B2SmyWldYa1fZvvKJ2aDw0AzazXi7O~oI6s7ayt3lWhCFAdaR4i6bACnaDag__"
tar -zxvf spaceranger-3.1.1.tar.gz

# Get the full path
pwd

# Expected output
# The path will change dpending on the compute setup you are using.
cp p -r /home/f06b22037/SSD2/utility/spaceranger /home/woodydrylab/FileShare
export PATH=/home/woodydrylab/FileShare/spaceranger/spaceranger-3.1.1:$PATH
export PATH=/home/f06b22037/SSD2/utility/spaceranger/spaceranger-3.1.1:$PATH

grep -c processor /proc/cpuinfo
grep MemTotal /proc/meminfo | cut -d ':' -f 2 | sed 's/^[ \t]*//'

spaceranger testrun --id=verify_install