#!/bin/bash      
#title			:ami_launcher.sh
#description	:This script will launch AMI to analyse RNA-seq from MMETSP.
#author			:Adam Monier
#date			:20140324
#version		:0.1    
#usage			:bash ami_launcher.sh strain_OR_transcriptome_ID config_file
#notes			:AWS CLI tools are required for this script.
#notes			:Modified from instance launcher from Shanrong Zhao's work (original project: stormbow) 
#bash_version	:GNU bash, version 3.2.51(1)-release (x86_64-apple-darwin13)
#===============================================================================

# Timestamp for temp and log files
now=$(date +"%Y-%m-%d")

# Command line arguments initialisation
data_id=$1
config_file=$2

########################################################
# Test for single or multiple transcriptomes 
########################################################
if [ "${data_id:0:6}" == "MMETSP" ]
then
	data_type="transcriptome"
	transcriptome_id=$1
	strain_id=`awk 'BEGIN{FS=","};$2=="'$transcriptome_id'"{print $1}' MMETSP_listing.csv`
	printf "$transcriptome_id" > ${data_id}_${now}.list	
elif [ "${data_id:0:6}" == "STRAIN" ]
then
	data_type="strain"
	strain_id=$1
	awk 'BEGIN{FS=","};$1=="'$strain_id'"{print $2}' MMETSP_listing.csv > ${data_id}_${now}.list
else 
	printf "Error: Data identifier must be either STRAIN or MMETSP\n"
	exit 1
fi

########################################################
# Get MMETSP information from MMETSP_listing.csv
########################################################
species=`awk 'BEGIN{FS=","};$1=="'$strain_id'"{print $6,$7,$8}' MMETSP_listing.csv|awk 'NR==1'|sed -e "s/\.//g" -e "s/ /_/g"`
phylum=`awk 'BEGIN{FS=","};$1=="'$strain_id'"{print $3}' MMETSP_listing.csv|awk 'NR==1'|sed -e "s/\.//g" -e "s/ /_/g"`
pinvest=`awk 'BEGIN{FS=","};$1=="'$strain_id'"{print $12}' MMETSP_listing.csv|awk 'NR==1'|sed -e "s/\.//g" -e "s/ /_/g"`
num_transcriptome=`awk 'END{print NR}' ${data_id}_${now}.list`

########################################################
# Get AWS information from ami configuration file (e.g., ami.cfg)
########################################################
SECURITY=`grep "^SECURITY=" $config_file | cut -f2 -d"="`
KEY=`grep "^KEY=" $config_file | cut -f2 -d"="`
REGION=`grep "^REGION=" $config_file | cut -f2 -d"="`
AK=`grep "^ACCESS_KEY=" $config_file | cut -f2 -d"="`
SK=`grep "^SECRET_KEY=" $config_file | cut -f2 -d"="`
INSTANCE=`grep "^INSTANCE=" $config_file | cut -f2 -d"="`
MAXPRICE=`grep "^MAXPRICE=" $config_file | cut -f2 -d"="`
S3DATA=`grep "S3DATA=" $config_file | cut -f2 -d"="`
# Hard-coded parameters for AWS AMI and S3 code repository (where RNAseq.sh is available)
AMI="ami-d58e70a2"
code_bucket="s3://mmetsp/code"
# Upload MMETSP list file to S3 code bucket
s3cmd --quiet -f put ${data_id}_${now}.list ${code_bucket}/${data_id}_${now}.list > /dev/null 2>&1

########################################################
# Prepare user-data-file to pass to instance
########################################################
printf "#!/bin/bash -ex\n" >> ${data_id}_${now}_data-user.sh

#test for instance type (t1 vs. m1) and mount ONE ephemeral drive
if [ "${INSTANCE:0:2}" == "m1" ]
then
	printf "umount /mnt\n" >> ${data_id}_${now}_data-user.sh
	printf "mount /dev/xvdb /mnt/data1\n\n" >> ${data_id}_${now}_data-user.sh
elif [ "${INSTANCE:0:2}" == "t1" ]
then
	printf "t1.micro instance; no ephemeral drive\n"
else 
	printf "Error: Only t1 and m1 instances are available\n"
	exit 1
fi

# Pass environment variables
printf "export INSTANCE=$INSTANCE \n" >> ${data_id}_${now}_data-user.sh
printf "export DATA_ID=$data_id \n" >> ${data_id}_${now}_data-user.sh
printf "export DATA_TYPE=$data_type \n" >> ${data_id}_${now}_data-user.sh
printf "export CODE_BUCKET=$code_bucket \n" >> ${data_id}_${now}_data-user.sh
printf "export SPECIES=$species \n" >> ${data_id}_${now}_data-user.sh
printf "export PHYLUM=$phylum \n" >> ${data_id}_${now}_data-user.sh
printf "export PINVEST=$pinvest \n\n" >> ${data_id}_${now}_data-user.sh

# Modify .s3cfg config file for s3cmd based on instance config file details	
printf "sed -i 's|access_key = .*|access_key = $AK|g' /home/ubuntu/.s3cfg\n" >> ${data_id}_${now}_data-user.sh
printf "sed -i 's|secret_key = .*|secret_key = $SK|g' /home/ubuntu/.s3cfg\n\n" >> ${data_id}_${now}_data-user.sh

# Retrieve RNAseq.sh wrapper from S3 code bucket	
printf "s3cmd --config=/home/ubuntu/\.s3cfg get ${code_bucket}/RNAseq\.sh  /home/ubuntu/RNAseq\.sh\n" >> ${data_id}_${now}_data-user.sh
printf "s3cmd --config=/home/ubuntu/\.s3cfg get ${code_bucket}/${data_id}_${now}\.list  /home/ubuntu/${data_id}_${now}\.list\n" >> ${data_id}_${now}_data-user.sh
printf "chmod +x /home/ubuntu/RNAseq\.sh\n" >> ${data_id}_${now}_data-user.sh
printf "bash /home/ubuntu/RNAseq.sh\n" >> ${data_id}_${now}_data-user.sh

########################################################
# Launch AMI using CLI tools
# Test for spot request
# If spot, wait for spot request fulfillment
########################################################
if [ "${MAXPRICE}" == "" ]
then
	# Launch non-spot instance
	ec2-run-instances $AMI --region $REGION --group $SECURITY --key $KEY --aws-access-key $AK --aws-secret-key $SK --instance-initiated-shutdown-behavior terminate -t $INSTANCE --user-data-file ${data_id}_${now}_data-user.sh -b "/dev/sdb=ephemeral0" > ${data_id}_${now}.instance_info
	inst_id=`grep "^INSTANCE" ${data_id}_${now}.instance_info | awk '{print $2}'`
else 
	# Request spot
	ec2-request-spot-instances $AMI -p $MAXPRICE -n 1 --region $REGION --group $SECURITY --key $KEY --aws-access-key $AK --aws-secret-key $SK --type one-time -t $INSTANCE --user-data-file ${data_id}_${now}_data-user.sh -b "/dev/sdb=ephemeral0" > ${data_id}_${now}.instance_info
	spot_id=`grep "^SPOTINSTANCEREQUEST" ${data_id}_${now}.instance_info | awk '{print $2}'`
	printf "Amazon Web Services spot $spot_id requested in $REGION region (based on ami: $AMI)\n"
	
	# Wait for spot request fulfillment
	inst_id=""
	wait_time=0
	while [ "${inst_id}" == "" ]
	do
		inst_id=`ec2-describe-spot-instance-requests $spot_id | grep "^SPOTINSTANCEREQUEST" | awk 'BEGIN{FS="\t"}{print $12}'`
		spot_status=`ec2-describe-spot-instance-requests $spot_id | grep "^SPOTINSTANCESTATUS" | awk 'BEGIN{FS="\t"}{print $2}'`
		printf "$wait_time seconds: Spot request status: $spot_status ...Waiting...\n"
		sleep 15
		wait_time=`expr $wait_time + 15`
	done
fi

########################################################
# Get instance info, edit spot/instance tags, print information
########################################################
ec2-create-tags $inst_id --tag "Launch time=$(date +'%Y-%m-%d-%S')" --tag "Name=$data_id" --tag "Phylum=$phylum" --tag "Species=$species" --tag "PI=$pinvest" --tag "MMETSP experiment(s)=$num_transcriptome" --aws-access-key $AK --aws-secret-key $SK > /dev/null
ec2-create-tags $spot_id --tag "Launch time=$(date +'%Y-%m-%d-%S')" --tag "Name=$data_id" --tag "Phylum=$phylum" --tag "Species=$species" --tag "PI=$pinvest" --tag "MMETSP experiment(s)=$num_transcriptome" --aws-access-key $AK --aws-secret-key $SK > /dev/null
printf "RNA-seq analysis for $species (strain ID: $strain_id; phylum: $phylum; PI: $pinvest), composed of $num_transcriptome transcriptome experiment(s)\n"
printf "Amazon Web Services instance $inst_id launched in $REGION region (based on ami: $AMI)\n"

# Clean up
yes|rm ${data_id}_${now}*

exit 1
