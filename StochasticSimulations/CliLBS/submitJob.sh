#!/bin/bash
lastLbs=$(ls -rt *lbs | tail -1)
number=$(echo ${lastLbs%.lbs}|cut -d'_' -f'2')
nextid=$(echo "${number}+1"| bc)
nextLbs=${lastLbs%${number}.lbs}$(printf %06d ${nextid}).lbs
echo $nextLbs
cp ${lastLbs} ${nextLbs}
sed -e s/#1/${nextLbs}/g queueTemplate.pbs > queue_${nextLbs%.lbs}.pbs
qsub queue_${nextLbs%.lbs}.pbs
