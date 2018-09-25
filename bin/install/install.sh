BUSCO_DATA="actinopterygii_odb9"

conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda env create --name snakemake_genome --file requirements.txt
source activate snakemake_genome

git clone https://gitlab.com/ezlab/busco.git ../busco
mkdir ../../data/busco
wget http://busco.ezlab.org/v2/datasets/$BUSCO_DATA.tar.gz -P ../../data/busco/
tar -xvf ../../data/busco/$BUSCO_DATA.tar.gz -C ../../data/busco/
rm ../../data/busco/$BUSCO_DATA.tar.gz

cp ../config.ini ../busco/config/
export AUGUSTUS_CONFIG_PATH="/home/chollenbeck/bin/pkgs/augustus-3.2.3/config/"
