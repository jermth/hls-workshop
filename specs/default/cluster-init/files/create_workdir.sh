#! /bin/bash

mkdir -p $HOME/genomics-workshop/
sourcedir="/avere/data/genomics_workshop_data"
if [ ! -L $HOME/genomics-workshop/bowtie2_index ];then
    ln -s $sourcedir/bowtie2_index $HOME/genomics-workshop/
fi

if [ ! -L $HOME/genomics-workshop/NA12787 ];then
    ln -s $sourcedir/NA12787 $HOME/genomics-workshop/
fi

cd $HOME/genomics-workshop/ 
for file in bowtie2.container.nf  bowtie2.local.nf  nextflow.config  nextflow.sge.config
do
    if [ ! -f $file ];then
        cp $sourcedir/$file $HOME/genomics-workshop/ 
    fi
done

cd $HOME
