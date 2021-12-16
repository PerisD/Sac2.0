#############################
# General options
#############################
genome_size = 72000000
genome = /home/GLBRCORG/pnavarro/GLBRC_SeptaHyb/analysis/assemblies/reference_genome/ParentalGenomes.fasta
out_dir = /home/GLBRCORG/pnavarro/GLBRC_SeptaHyb/analysis/assemblies/STRAINAME	        # all outputs will be written to the folder "./iWGS_test"
threads = 4				        	# number of CPUs to use (default: 1)
memory = 32				            	# number of GBs of memory to use (default: 8)

#############################
# Library options
#############################
library = STRAINAME,PE,5,150,280,50				



#############################
# Simulator options
#############################
pIRS.error_rate = 				        # substitution error rate: 0, 1, or 0.0001-0.63 (default: 1 - indicate that the default setting of pIRS should be used)
pIRS.error_profile = 		        		# the base-calling profile for simulating substitution-error and quality score (default: the default profile of pIRS)
pIRS.gc = 			        		# whether to simulate GC bias: 1 - yes; 0 - no (default: 1)
pIRS.gc_profile = 			        	# the GC content-coverage file for simulating GC bias (default: the default profile of pIRS)
pIRS.indel = 				        	# whether to simulate indel errors: 1 - yes; 0 - no (default: 1)
pIRS.indel_profile =			 	    	# the InDel-error profile for simulating InDel-error (default: the default profile of pIRS)

ART.qual_shift1 = 		        		# the amount to shift the quality score of all first-reads (default: 0)
ART.qual_shift2 = 				        # the amount to shift the quality score of all second-reads (default: 0)
ART.qual_profile1 = 		      			# the quality profile of first-reads (default: the default profile of ART)
ART.qual_profile2 = 		        		# the quality profile of second-reads (default: the default profile of ART)
ART.ins_rate1 = 		        		# the insertion rate of first-reads (default: 0.00009)
ART.ins_rate2 = 		        		# the insertion rate of second-reads (default: 0.00015)
ART.del_rate1 = 		        		# the deletion rate of first-reads (default: 0.00011)
ART.del_rate2 = 		        		# the deletion rate of second-reads (default: 0.00023)

PBSIM.model_qc = /home/GLBRCORG/pnavarro/software/iWGS_v9/tools/PBSIM_model/model_qc_clr-alyrata				        # the model of quality code for simulating read accuracy (default: the default profile of PBSIM)
PBSIM.ratio = 					        # the ratio of substitution:insertion:deletion errors (default: 10:60:30)
PBSIM.accuracy_max = 					# the maximum read accuracy (default: 0.90)
PBSIM.accuracy_min = 				        # the minimum read accuracy (default: 0.75)
PBSIM.length_mean = 			        	# the mean of read length (default: 3000)
PBSIM.length_sd = 			        	# the standard deviation of read length (default: 2300)
PBSIM.length_max = 		        		# the maximum read length (default: 25000)
PBSIM.length_min = 					# the minimum read length (default: 100)

#############################
# Quality control options
#############################
QC = STRAINAME						# whether to perform QC on selected libraries (provide a list of library names separated by comma) or all libraries ("all")

Trimmomatic.trailing = 3				# the quality score cutoff for Trimmomatic quaulity-based trimming from 3' end of reads
Trimmomatic.adapters = 		# the file containing the list of adapter sequence files used for Trimmomatic adapter trimming 				
Trimmomatic.minlen = 25					# the minimum read length after Trimmomatic trimming

NextClip.adapter = 					# the adapter sequence for NextClip adapter trimming
NextClip.minlen = 25					# the minimum read length after NextClip trimming

Correction.tool = Lighter				# the error correction tool to be used (default: "Lighter"; also support "Quake")
Correction.kmer = 					# the k-mer size used for error correction

#############################
# Assembly protocol options
#############################
protocol = STRAINAMESO,SOAPdenovo2,STRAINAME
protocol = STRAINAMEVe,Velvet,STRAINAME
protocol = STRAINAMESP,SPADES,STRAINAME
protocol = STRAINAMEAB,ABYSS,STRAINAME
protocol = STRAINAMESG,SGA,STRAINAME
protocol = STRAINAMEMa,MaSuRCA,STRAINAME
protocol = STRAINAMEDI,DISCOVAR,STRAINAME
protocol = STRAINAMEPl,Platanus,STRAINAME
protocol = STRAINAMEDip,dipspades,STRAINAME


#############################
# Assembler options
#############################
ABYSS.kmer = 				        	# default: 0 - to be estimated using KmerGenie 
ABYSS.option = "l=1 n=5 s=100"				# GAGE-B recipe

ALLPATHS.ploidy = 2				        # default: 1; set to 2 for heterozygous assembly

CA.pbCNS =					        # default: 1;
CA.sensitive =					       	# default: 0; set to 1 for lower-quality PacBio data

Canu.sensitive =					# default: 0; set to 1 for lower-quality PacBio data

MaSuRCA.kmer =			        		# default: 0 - to be determined automatically by MaSuRCA

Meraculous.kmer = 			        	# default: 0 - to be determined automatically by Meraculous
Meraculous.diploid = 				    	# default: 0; set to 1 for diploid genome sequencing data

SPAdes.kmer = 				        	# default: 0 - to be estimated using KmerGenie, and used in addition to "multi-kmer" (if turned on)
SPAdes.multi-kmer = 				    	# default: 1 - to use multiple k-mer sizes
SPAdes.option = "--only-assembler"			# set it empty " " to enable the error-correction module

SOAPdenovo2.kmer = 			        	# default: 0 - to be estimated using KmerGenie
SOAPdenovo2.option = "-F -R -E -w -u"			# GAGE-B recipe

Velvet.kmer = 				        	# default: 0 - to be estimated using KmerGenie
Velvet.option = "-exp_cov auto -scaffolding yes"	# GAGE-B recipe

#############################
# Evaluation options
#############################
QUAST.eukaryote = 1	    				# whether the reference genome is eukaryotic: 1 - yes; 0 - no (default: 1)
QUAST.gage = 1			        		# whether to generate GAGE report: 1 - yes; 0 - no (default: 1)
QUAST.gene = 				        	# gene annotations to be used for evaluation (default: NA)

REAPR.libs = 				                # the libraries to be used for REAPR evalution

#############################
# Executable options
#############################
bin.ART = /home/GLBRCORG/pnavarro/software/ART/art_illumina  ############
bin.pIRS = /opt/bifxapps/pirs/pirs
bin.PBSIM = /opt/bifxapps/pbsim-1.0.3/Linux-amd64/bin/pbsim
bin.Trimmomatic = /opt/bifxapps/trimmomatic/trimmomatic-0.30.jar
bin.NextClip = /home/GLBRCORG/pnavarro/software/nextclip-master/bin/nextclip
bin.Lighter = /home/GLBRCORG/pnavarro/software/Lighter-master/lighter
bin.Quake = /home/GLBRCORG/pnavarro/Quake/bin/quake.py
bin.KmerGenie = /opt/bifxapps/kmergenie-1.6950/kmergenie
bin.ABYSS = /opt/bifxapps/abyss-1.5.2/bin/abyss-pe
bin.ALLPATHS = /opt/bifxapps/allpathslg/bin/RunAllPathsLG
bin.BLASR = /home/GLBRCORG/pnavarro/software/blasr
bin.PBcR = /mnt/bigdata/processed_data/pnavarro/iwgs_data/wgs-8.3rc1/Linux-amd64/bin/PBcR
bin.PBDAGCON = /opt/bifxapps/pbdagcon/pbdagcon
bin.runCA = /mnt/bigdata/processed_data/pnavarro/iwgs_data/wgs-8.3rc1/Linux-amd64/bin/runCA
bin.DISCOVAR = /opt/bifxapps/discovarexp/bin/DiscovarExp
bin.MaSuRCA = /opt/bifxapps/masurca/bin/masurca
#bin.Meraculous = /home/GLBRCORG/jacek.kominek/software/meraculous/bin/run_meraculous.sh #Jacek comments it is not improving any assembly
bin.Platanus = /opt/bifxapps/bin/platanus
bin.SGA = /opt/bifxapps/sga/bin/sga
bin.SOAPdenovo2 = /opt/bifxapps/SOAPdenovo2/SOAPdenovo-127mer
bin.SPAdes = /opt/bifxapps/SPAdes-3.5.0/bin/spades.py
bin.dipSPAdes = /opt/bifxapps/SPAdes-3.5.0/bin/dipspades.py
bin.velveth = /opt/bifxapps/bin/velveth
bin.velvetg = /opt/bifxapps/bin/velvetg
bin.QUAST = /home/GLBRCORG/pnavarro/software/quast-3.2/quast.py
bin.REAPR = /home/GLBRCORG/jacek.kominek/software/iWGS_v1.0/tools/REAPR/reapr
bin.bank-transact = /opt/bifxapps/amos/bin/bank-transact
bin.BWA = /opt/bifxapps/bin/bwa
bin.SAMtools = /home/GLBRCORG/jacek.kominek/software/iWGS_v1.0/tools/dependencies/samtools
bin.FASTX = /home/GLBRCORG/jacek.kominek/software/iWGS_v1.0/tools/dependencies/fastx_reverse_complement
bin.Metassembler = /home/GLBRCORG/pnavarro/software/Metassembler/bin/metassemble
bin.FALCON = /home/GLBRCORG/pnavarro/software/FALCON/fc_env/bin/fc_run.py
bin.Canu = /home/GLBRCORG/pnavarro/software/canu/bin/canu
#bin.DBG2OLC = /home/GLBRCORG/pnavarro/software/
#bin.Sparc = /home/GLBRCORG/pnavarro/software/
bin.Minia = /opt/bifxapps/minia/minia