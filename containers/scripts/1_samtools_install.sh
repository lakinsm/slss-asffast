#!/usr/bin/env bash
# written by oleg osipenko
# 10012020
#

# ----- array of needed apps -----
packages[0]='build-essential'
packages[1]='git'
packages[2]='libz-dev'
packages[3]='libbz2-dev'
packages[4]='zlib1g-dev'
packages[5]='libncurses5-dev'
packages[6]='libncursesw5-dev'
packages[7]='liblzma-dev'
packages[8]='gfortran'
packages[9]='curl'
packages[10]='wget'

# ----- function to check if package exists -----
function package_exists(){
	dpkg -s "$1" &> /dev/null
	return $?
}

# ----- check if package exists for each element in array -----
for pkg in "${packages[@]}"
do
	if ! package_exists $pkg ; then
		echo "$pkg is not isntalled, installing now" 
		apt-get install -y $pkg
	else
		echo "$pkg is installed ... good to go"
	fi
done

# -----  set directory for biotools and start cloning -----
[ ! -d /opt/biotools ]
	mkdir /opt/biotools
cd /opt/biotools || exit

wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 
tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9 || exit
make
cd ../

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9 || exit
make
cd ../

wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9 || exit
make
cd ../

git clone https://github.com/lh3/minimap2.git
cd minimap2 || exit
make
cd ../

# -----  rename for congruency with python script -----
# mv ncbi* MCBI
mv ncbi-blast-2.10.1+ NCBI 
mv bcftools-1.9 bcftools
mv htslib-1.9 htslib
mv samtools-1.9 samtools
cp bcftools/misc/vcfutils.pl /opt/biotools/
cd /usr/bin
ln -s /opt/biotools/minimap2/minimap2
ln -s /opt/biotools/bcftools/bcftools
ln -s /opt/biotools/vcfutils.pl
ln -s /opt/biotools/samtools/samtools
