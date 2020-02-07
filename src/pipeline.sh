set -e
type='EBV1'
referenceEBV='/home/cata/EBV/Genoma_referencia/'$type'/*.fa'
reference='/home/cata/EBV/Genoma_referencia/'$type'/*NNN.fasta'
regionfile='/home/cata/EBV/coordenadas/'$type'/Coordenadas.'$type

#samples=($(ls /home/cata/Escritorio/13-09-2019/Secuencias_crudas/))
#declare -a samples=("M1" "M2" "M3" "M4" "M5" "M6" "M7" "M8" "M9" "M10" "M11" "M13" "M14" "M15" "M16" "M12")
#declare -a samples=("M1")
samples='M1'
outpath='/home/cata/EBV/20190913/'$type
outtrimmed='/home/cata/EBV/20190913'
#sudo mkdir $outpath -p

for s in "${samples[@]}"
	do
	echo $s
	#mkdir $outpath/$s -p
	#mkdir $outtrimmed/$s -p
	#zcat /home/cata/Escritorio/13-09-2019/Secuencias_crudas/$s/*/*/$s*"R1"* |gzip > $outtrimmed/$s/$s.R1.fq.gz
	#zcat /home/cata/Escritorio/13-09-2019/Secuencias_crudas/$s/*/*/$s*"R2"* |gzip > $outtrimmed/$s/$s.R2.fq.gz

	#fastp --disable_adapter_trimming -i $outtrimmed/$s/$s.R1.fq.gz -I $outtrimmed/$s/$s.R2.fq.gz -o $outtrimmed/$s/$s.R1.trimmed.fq.gz \
	#-O $outtrimmed/$s/$s.R2.trimmed.fq.gz -h \
	#$outtrimmed/$s/$s.report_fastp.html -j $outtrimmed/$s/$s.report_fastp.json \
	#--trim_front1=15 --trim_tail1=5 -e 28 --length_required 100;

	## map to reference
	#bwa mem -K 100000000 -v 1 -t 4 $referenceEBV <(zcat $outtrimmed/$s/$s.R1.trimmed.fq.gz) <(zcat $outtrimmed/$s/$s.R2.trimmed.fq.gz) | samtools view -b - > $outpath/$s/$s.bam
	
	## remove duplicates and unmapped reads
	#samtools view -b -F 1548 $outpath/$s/$s.bam > $outpath/$s/$s.mapped.bam  
	
	## sort bam
	#samtools sort $outpath/$s/$s.mapped.bam > $outpath/$s/$s.mapped.sorted.bam
	
	# drop out repetitive regions (according to coords file /home/cata/). 
	#samtools view -bL $regionfile $outpath/$s/$s.mapped.sorted.bam > $outpath/$s/$s.mapped.sorted.withoutrep.bam

	## create bam index
	#samtools index $outpath/$s/$s.mapped.sorted.withoutrep.bam
	
	# compute statistics for whole bam
	#samtools stats $outpath/$s/$s.bam > $outpath/$s/$s.stats

	# compute statistics for final bam (without repetitive regions nither unmapped reads)
	#samtools stats $outpath/$s/$s.mapped.sorted.withoutrep.bam > $outpath/$s/$s.mapped.sorted.withoutrep.stats
	
	## get vcf
	#bcftools mpileup -f $referenceEBV $outpath/$s/$s.mapped.sorted.withoutrep.bam |bcftools call -mv --ploidy 1 -o $outpath/$s/$s.calls.vcf
	###bcftools index $outpath/$s/$s.calls.vcf.gz
	###gzip -dk $outpath/$s/$s.calls.vcf.gz
	
	#bgzip -c $outpath/$s/$s.calls.vcf > $outpath/$s/$s.calls.vcf.gz
	#tabix $outpath/$s/$s.calls.vcf.gz


	#Cat VCF
	interval='/home/cata/EBV/coordenadas/EBV1/Coordenadas.EBV1.bed' ## hay que generarlo from scrach a partir de 
	zgrep '^#' $outpath/$s/$s.calls.vcf.gz > header
	intersectBed -wa -a $outpath/$s/$s.calls.vcf.gz -b $interval > body.vcf  ### 
	cat header body.vcf > $outpath/$s/$s.NonRep.calls.vcf
	
	bgzip -c $outpath/$s/$s.NonRep.calls.vcf > $outpath/$s/$s.NonRep.calls.vcf.gz
	tabix $outpath/$s/$s.NonRep.calls.vcf.gz
	rm header
	rm body.vcf

	#Coverange_0
	bedtools genomecov -ibam $outpath/$s/$s.mapped.sorted.bam -bga | awk '$4==0'| less > $outpath/$s/$s.Cob0.bed

	#Deletion
	vcf2bed -n < $outpath/$s/$s.NonRep.calls.vcf > $outpath/$s/$s.deletion.bed
	less $outpath/$s/$s.deletion.bed | awk '{print $1, $2, $3}' > $outpath/$s/$s.deletion.c.bed
	less $outpath/$s/$s.deletion.c.bed |tr ' ' '\t' > $outpath/$s/$s.deletion.bed
	rm $outpath/$s/$s.deletion.c.bed

	
	## get consensus sequence
	bedtools subtract -a $outpath/$s/$s.Cob0.bed -b $outpath/$s/$s.deletion.bed > $outpath/$s/$s.COB-DEL.bed
	modified_reference=$outpath/$s/$s'_modified_reference.fa'
	/home/cata/EBV/src/mask_reference.py -r $reference -c $outpath/$s/$s.COB-DEL.bed -o $modified_reference
	bcftools consensus -f $modified_reference $outpath/$s/$s.NonRep.calls.vcf.gz > $outpath/$s/$s.nonrep.consensus.fa
	
done
 
