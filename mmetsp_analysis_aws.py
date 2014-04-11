#!/usr/bin/python      
#title			:mmetsp_analysis.py
#description	:This script will launch a de novo MMETSP transcriptome assembly
#author			:Adam Monier (marinebugs~at~gmail~dot~com)
#date			:20140408
#version		:0.2    
#usage			:python mmetsp_analysis.py -i [MMETSP identifier] > log 2>&1 &
#notes			:AWS CLI tools are required for this script if backup in the cloud.
#notes			:Bioinformatics 3rd party tools: trimmomatic-0.32, Trinity, deconseq
#notes			:blast+, fastQC, prinseq
#notes			:All requirements fulfilled on AWS AMI: xxxxxx
#python_version	:Python 2.7.6, Anaconda 1.9.1 (x86_64)
#===============================================================================

# symbolic links: deconseq.pl, Trinity.pl, interleave, split, run_RSEM_align_n_estimate.pl, pick_isoform_trinity_RSEM.py
import os, sys, glob, subprocess, getopt
from ftplib import FTP


APPS_DIR = '/home/ubuntu/apps'
ADAPT = '/home/ubuntu/apps/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa'
#AWS_EPH = '/mnt/data1' # instance ephemeral mount
TMP_DIR = '/home/ubuntu/tmp'
DECONSEQ_DB = 'plast'
BLAST_CHIM_DB = '/home/ubuntu/data/protist/protistDB.pep'
HEADER = '\n'+120*'#'+'\n'+25*' '

##### Get CPU and RAM values
cmd = "awk 'NR==1{print $2/1048576-0.5}' /proc/meminfo"  # keep some mem.
mem_proc = subprocess.check_output(cmd, shell=True)
MEM_TOT = mem_proc[0]
cpu_proc = subprocess.Popen('nproc', stdout=subprocess.PIPE) # assign CPU number from nproc
CPU_TOT = cpu_proc.stdout.read().rstrip('\n')
#####

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def buildProjectDir(m_id):
	"""docstring"""
	print HEADER + 'Building project directory structure for ' + m_id + HEADER
	subprocess.call('mkdir ' + TMP_DIR + '/' + m_id, shell=True)
	dirs1 = ['reads', 'stats', 'assemblies', 'deconseq', 'proteins', 'annotation']
	for dir1 in dirs1:
		subprocess.call('mkdir ' + TMP_DIR + '/' + m_id + '/' + dir1, shell=True)
	dirs2 = ['ncgr','trim','trim.diginorm','trim.diginorm.deconseq','trinity','trinity.pickH',\
	'trinity.pickH.reduced','trinity.pickH.reduced.chimera']
	for dir2 in dirs2:
		subprocess.call('mkdir ' + TMP_DIR + '/' + m_id + '/stats/' + dir2, shell=True)
	
	
def downloadData(m_id,type):
	""" Retrieves MMETSP data from CAMERA ftp. For a listing of available transcriptomes \
	check http://bit.ly/1e2ZqvC """
	ftp = FTP('portal.camera.calit2.net','anonymous','anonymous@sunet.se')
	ftpdir = 'ftp-links/cam_datasets/projects/additional/CAM_P_0001000'	
	if type == 'reads':
		ftpfile = ftpdir + '/read/{}.fastq.tar'.format(m_id)
		ftpout = m_id + '.ncgr.fq.tar'	
	elif type == 'annotation':
		ftpfile = ftpdir + '/annotation/{}.annot.tgz'.format(m_id)
		ftpout = m_id + '.ncgr.annot.tgz'
	elif type == 'peptides':
		ftpfile = ftpdir + '/assemblies/{}.pep.fa.gz'.format(m_id)
		ftpout = m_id + '.ncgr.pep.fa.gz'
	elif type == 'contigs':
		ftpfile = ftpdir + '/assemblies/{}.nt.fa.gz'.format(m_id)
		ftpout = m_id + '.ncgr.nt.fa.gz'	
	elif type == 'cds':
		ftpfile = ftpdir + '/assemblies/{}.cds.fa.gz'.format(m_id)
		ftpout = m_id + '.ncgr.cds.fa.gz'
	else:
		print 'wrong data type. Valid types: reads, annotation, peptides, contigs, cds'
	print HEADER + 'Downloading' + ftpfile + HEADER
	ftp.retrbinary('RETR ' + ftpfile, open(ftpout, 'wb').write)
	ftp.quit()
	


def callTrimmomatic(pe1,pe2,trim_pe1,trim_pe2,trim_se1,trim_se2,adapters,qtrail,headclip):
	""" Launches Trimmomatic v0.32 (http://bit.ly/1hvbUqH) on PE fastq files.
	Detects and removes adapters from given fasta file (variable: ADAPT),
	removes first 13-bp (see paper http://1.usa.gov/1hQpYzm)"""
	print HEADER +'Trimmomatic on ' + pe1 + ' and ' + pe2 + HEADER
	trimdir = APPS_DIR + '/Trimmomatic-0.32'
	subprocess.call('java -jar {0}/trimmomatic-0.32.jar PE {1} {2} {3} {4} {5} {6} \
	ILLUMINACLIP:{7}:2:30:10 TRAILING:{8} HEADCROP:{9}'.format(trimdir,pe1,pe2,trim_pe1,\
	trim_se1,trim_pe2,trim_se2,adapters,qtrail,headclip), shell=True)
	
def callDigiNorm(pe,Cval,kval,Nval,xval):
	""" single-pass digital normalization (arXiv: http://bit.ly/1ifqPGs) using recommended \
	parameters.
	Script normalize-by-median.py is part of the khmer suite (GitHub: http://bit.ly/PLHMB8)."""
	print HEADER +'Digital normalization on ' + pe + HEADER
	subprocess.call('normalize-by-median.py -C {0} -k {1} -N {2} -x {3} {4}'\
	.format(Cval,kval,Nval,xval,pe), shell=True)
	 # might need to add APPS_DIR + bin/
	

def callDeconseq(pe,id,dbs,outdir):
	"""docstring"""
	print HEADER +'Deconseq on ' + pe + HEADER
	subprocess.call('deconseq.pl -f {0} -id {1}.trim.pe.diginorm.deconseq -keep_tmp_files\
	-dbs {2} -out_dir {3} 2>/dev/null'.format(pe,id,dbs,outdir), shell=True)
	
def fancyPlots(pe,outdir):
	"""docstring"""
	print 120*'#'+'\n'+25*' '+'Plots for ' + pe + '\n' + 120*'#'
	plotlist = ['plot_scores','plot_nucleotide_distribution','gc_distribution']
	for plot in plotlist:
		if plot in ('plot_scores','plot_distribution','plot_nucleotide_distribution'):
			subprocess.call("read_fastq -i {0} | {1} -T '{0} {1}' -t svg -o {2}/{0}.{1}.svg \
			-x".format(pe,plot,outdir), shell=True)
		elif plot == "gc_distribution":
			subprocess.call("read_fastq -i {0} | analyze_gc | bin_vals -b 5 -k GC% | \
			plot_distribution -k GC%_BIN -T '{0} {1}' -t svg -o {2}/{0}.gc_plot.svg -x"\
			.format(pe,plot,outdir), shell=True)
		elif plot == "plot_distribution":
			subprocess.call("read_fastq -i {0} | {1} -k SEQ_LEN -T '{0} {1}' -t svg -o {2}/{0}.{1}.svg \
			-x".format(pe,plot,outdir), shell=True)
		else:
			print "Error in biopiece plot type"
			
def callFastQC(pe,outdir):
	"""docstring"""
	print HEADER +'FastQC analysis for ' + pe + HEADER
	subprocess.call('fastqc {0} -o {1}'.format(pe,outdir), shell=True)
	
	
def callTrinityAssembler(pe1,pe2,mem,cpu,out):
	"""docstring"""
	print HEADER + 'Trinity assembly using ' + pe1 +' ' + pe2 + HEADER
	subprocess.call('ulimit -s unlimited', shell=True)
	subprocess.call('Trinity.pl --full_cleanup --seqType fq --SS_lib_type RF --left {0} --right {1} \
	--JM {2} --CPU {3} --output {4}'.format(pe1,pe2,mem,cpu,out), shell=True)

	
def assemblyReduction(id,contig,petrim1,petrim2,cpu,outdir):
	#subprocess.call('run_RSEM_align_n_estimate.pl --transcripts /home/ubuntu/tmp/dummy/assemblies/dummy.trinity.contig.nt.fa --seqType fq --left /home/ubuntu/tmp/dummy/reads/dummy.trim.pe1.fq --right /home/ubuntu/tmp/dummy/reads/dummy.trim.pe2.fq --thread_count 2 --output_dir /home/ubuntu/tmp/dummy/assemblies', shell=True)
	cw = os.getcwd()
	print HEADER + 'RSEM using ' + petrim1 +' and ' + petrim2 + ' on ' + contig + HEADER
	subprocess.call('run_RSEM_align_n_estimate.pl --transcripts {0} --seqType fq --left {1} --right {2} --thread_count {3} --output_dir {4}'.format(contig,petrim1,petrim2,cpu,outdir), shell=True)
	print HEADER + 'Picking isoforms for ' + contig + HEADER
	subprocess.call('pick_isoform_trinity_RSEM.py {0}/RSEM.isoforms.results {1}'.format(outdir,contig), shell=True)
	os.rename('{0}/{1}.trinity.contig.nt.fa.exemplar'.format(outdir,id),'{0}/{1}.trinity.pickH.contig.nt.fa'.format(outdir,id))
	print HEADER + 'CAP3 on ' + '{0}/{1}.trinity.pickH.contig.nt.fa'.format(outdir,id) + HEADER
	subprocess.call('cap3 {0}/{1}.trinity.pickH.contig.nt.fa -o 200 -p 99'.format(outdir,id), shell=True)
	os.rename('{0}/{1}.trinity.pickH.contig.nt.fa.cap.contigs'.format(outdir,id),'{0}/{1}.trinity.pickH.cap3.contig.nt.fa'.format(outdir,id))
	print HEADER + 'CD-HIT-EST on ' + '{0}/{1}.trinity.pickH.cap3.contig.nt.fa'.format(outdir,id) + HEADER
	subprocess.call('cd-hit-est -i {0}/{1}.trinity.pickH.cap3.contig.nt.fa -o {0}/{1}.trinity.pickH.cap3.cdhit.contig.nt.fa -c 0.99'.format(outdir,id),shell=True)
	#subprocess.call('pick_isoform_trinity_RSEM.py {0}/RSEM.isoforms.results {1}'.format(outdir,contig))

def chimeraCheck(id,query,output,outdir,db,maxseq,eval,minsize):
	"""docstring"""
	print HEADER +'BLASTX of ' + query + ' vs. ' + db + HEADER
	outformat = "'6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'"
	print 'blastx -db {0} -query {1} -evalue {2} -outfmt {3} -out {4} -num_threads={5} -max_target_seqs {6}'.format(db,query,eval,outformat,output,CPU_TOT,maxseq)
	subprocess.call('blastx -db {0} -query {1} -evalue {2} -outfmt {3} -out {4} -num_threads={5} -max_target_seqs {6}'.format(db,query,eval,outformat,output,CPU_TOT,maxseq), shell=True)
	print 120*'#'+'\n'+25*' '+'Detect chimera from ' + output + '\n' + 120*'#'
	subprocess.call('detect_chimera_from_blastx.py ' + outdir, shell=True)
	subprocess.call('cp assemblies/{0}.trinity.pickH.cap3.cdhit.contig.blastx.cut assemblies/{0}.blastx.cut'.format(id),shell=True) # copy and rename cut output for next step
	subprocess.call('cp {0} {0}.before_cut'.format(query),shell=True)
	print 120*'#'+'\n'+25*' '+'Cut chimeras from ' + query  + '\n' + 120*'#'
	subprocess.call('cut_chimera_from_blastx.py {0} {0} {1}'.format(outdir,minsize),shell=True)
	
# def assemblyStat():
	#"""docstring"""
	
# transdecoder

# buildResultLog()

def cleanUp(id):
	"""docstring"""
	subprocess.call('rm stats/*/*.zip',shell=True)
	subprocess.call('rm -r assemblies/RSEM.*',shell=True)
	subprocess.call('rm assemblies/*.fa.*',shell=True)
	subprocess.call('rm assemblies/{0}.blastx.cut'.format(id),shell=True)
	
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:",["mmetsp="])
	except getopt.GetoptError:
		print 'mmetsp_analysis.py -i <mmetsp identifier>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'mmetsp_analysis.py -i <mmetsp identifier>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			M_ID = arg.split('_')[0] # need to remove potential "_2" from mmetsp ID

	print M_ID
	buildProjectDir(M_ID) #will need to add AWS_EPH/
	with cd(TMP_DIR + '/' + M_ID): # project context
		#downloadData(M_ID,'reads') # download read tar archive
		subprocess.call('cp /home/ubuntu/apps/MMETSP/dummy.fastq.tar {0}/{1}'\
		.format(TMP_DIR,M_ID), shell=True) # using dummy instead of downloaded MMETSP data
		
		subprocess.call('tar xvf {0}.fastq.tar; rm {0}.fastq.tar'.format(M_ID), \
		shell=True) # untar read archive
		
		for filename in glob.glob('*.fastq.gz'): # rename fq files and gunzip
			if filename[-10] == '1':
				os.rename(filename,'reads/' + M_ID + '.ncgr.pe1.fq.gz') # rename file
				subprocess.call('gunzip -f reads/' + M_ID + '.ncgr.pe1.fq.gz', shell=True) # unzip fq
				ncgr_fq1 = M_ID + '.ncgr.pe1.fq' # assign vairable name to ncgr fq file
			elif filename[-10] == '2':
				os.rename(filename,'reads/' + M_ID + '.ncgr.pe2.fq.gz')
				subprocess.call('gunzip -f reads/' + M_ID + '.ncgr.pe2.fq.gz', shell=True) 
				ncgr_fq2 = M_ID + '.ncgr.pe2.fq'
			else: # in case fq is absent or bad name
				print 'Error in main(): bad file names in renaiming NCGR fq files'
				
		with cd(TMP_DIR + '/' + M_ID + '/reads'): # reads context
			callTrimmomatic(ncgr_fq1,ncgr_fq2,M_ID + '.trim.pe1.fq',M_ID + '.trim.pe2.fq',M_ID +\
			 '.trim.se1.fq',M_ID + '.trim.se2.fq', ADAPT, 30, 13) # adapter+qtrim+headclip
			
			subprocess.call('interleave-reads.py -o {0}.trim.pe.fq {0}.trim.pe1.fq {0}.trim.pe2.fq'\
			.format(M_ID), shell=True) # Interleave prior to DigitNorm, dump SE
			
			subprocess.call('interleave-reads.py -o {0}.ncgr.pe.fq {0}.ncgr.pe1.fq {0}.ncgr.pe2.fq'\
			.format(M_ID), shell=True) # Interleave ncgr reads because of subsequent fastQC
			
			callDigiNorm(M_ID + '.trim.pe.fq','20','20','4',str(float(MEM_TOT)/4) + 'e9') # DigiNorm on clean PE single fq file
			os.rename(M_ID + '.trim.pe.fq.keep', M_ID + '.trim.diginorm.pe.fq') # rename DigiNorm fq output
			ncgr_fq_dn = 'reads/{}.trim.diginorm.pe.fq'.format(M_ID)
			
		callDeconseq(ncgr_fq_dn,M_ID,DECONSEQ_DB,'deconseq/') # Deconseq on DigiNorm clean fq
		os.rename('deconseq/' + M_ID + '.trim.pe.diginorm.deconseq_clean.fq','reads/' \
		+ M_ID + '.trim.diginorm.deconseq.pe.fq') # move deconseq fq output to reads folder
		os.rename('deconseq/' + M_ID + '.trim.pe.diginorm.deconseq_cont.fq','deconseq/' \
		+ M_ID + '.trim.diginorm.deconseq_contaminant.fq') # rename deconseq contaminant fq
			
		#compute read stats with loop list (ncgr, trim, diginorm, deconseq)
		# need to interleave ncgr reads
		with cd(TMP_DIR + '/' + M_ID + '/reads'):
			readlist= ['ncgr','trim','trim.diginorm','trim.diginorm.deconseq']
			for readtype in readlist:
				callFastQC('{}.{}.pe.fq'.format(M_ID,readtype),'../stats/' + readtype)
				#fancyPlots('{}.{}.pe.fq'.format(M_ID,readtype),'../stats/' + readtype)
			subprocess.call('split-paired-reads.py '+ M_ID + '.trim.diginorm.deconseq.pe.fq',\
			shell=True) # split paired fq before Trinity, then rename files:
			os.rename(M_ID + '.trim.diginorm.deconseq.pe.fq.1',M_ID + '.trim.diginorm.deconseq.pe1.fq')
			os.rename(M_ID + '.trim.diginorm.deconseq.pe.fq.2',M_ID + '.trim.diginorm.deconseq.pe2.fq')
			callTrinityAssembler(M_ID + '.trim.diginorm.deconseq.pe1.fq',M_ID + '.trim.diginorm.deconseq.pe2.fq',\
			MEM_TOT + 'G',CPU_TOT,TMP_DIR + '/' + M_ID + '/assemblies/' + M_ID)	
			os.rename('../assemblies/' + M_ID + '.Trinity.fasta','../assemblies/' + M_ID + '.trinity.contig.nt.fa') 			
		with cd(TMP_DIR + '/' + M_ID):
			
			assemblyReduction(M_ID,'assemblies/' + M_ID + '.trinity.contig.nt.fa','reads/' + M_ID + '.trim.pe1.fq','reads/' + M_ID + '.trim.pe2.fq',CPU_TOT,'assemblies')
			subprocess.call('cp assemblies/{0}.trinity.contig.nt.fa assemblies/{0}.trinity.pickH.cap3.cdhit.contig.nt.fa'.format(M_ID), shell=True) #TEST
			chimeraCheck(M_ID,'assemblies/'+ M_ID + '.trinity.pickH.cap3.cdhit.contig.nt.fa','assemblies/'+ M_ID + '.trinity.pickH.cap3.cdhit.contig.blastx','assemblies',BLAST_CHIM_DB,100,0.01,1)
			
			
			
			cleanUp(M_ID)
			
			
			#download remaining data with a list loop
			
			# after final assembly rename contig ids with prinseq (starting w/ M_ID)
	
if __name__ == "__main__":
	main(sys.argv[1:])

