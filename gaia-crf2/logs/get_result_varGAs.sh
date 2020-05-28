#!/bin/bash

cat=('GA3.0' 'GA3.5' 'GA4.0' 'GA4.5' 'GA5.0' 'GA5.5' 'GA6.0' 'GA6.5' 'GA7.0' 'GA7.5' 'GA8.0' 'GA8.5' 'GA9.0' 'GA9.5')

for cati in ${cat[@]}
do
 	grep "used_for_grep_all" ${cati}_15.5_GaiaDR2_vsh_param.log >> vsh01_varGAs.log
 	grep "used_for_grep_GA" ${cati}_15.5_GaiaDR2_vsh_param.log >> GA_varGAs.log
done
