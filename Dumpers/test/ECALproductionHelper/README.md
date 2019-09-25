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

