# Examples for running
Productions conventions  https://docs.google.com/spreadsheets/d/1rJ4FppHxdD6TQpCHnnnX0QD0Y0d4aTEN6ibgTVn7Vyg/edit#gid=0

## EB no PU
python prodHelper.py -v V01_v01 -n 15000 -g closeEcal -d EB --pu noPU --domedium --domultithread 

## EB w/ PU
python prodHelper.py -v V01_v03 -n 15000 -g closeEcal -d EB --pu wPU --dolong --domultithread 

## EE no PU
python prodHelper.py -v V01_v02 -n 30000 -g closeEcal -d EE --pu noPU --dolong --domultithread

## EE w/ PU
python prodHelper.py -v V01_v04 -n 30000 -g closeEcal -d EE --pu wPU --dolong --domultithread

# Example for post-production
Production label as per conventions above

python postProdHelper.py --pl photon_Et1to100GeV_closeEcal_EE_noPU_pfrh1.0_seed3.0_V01_v52_n15000
