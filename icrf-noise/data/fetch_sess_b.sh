#########################################################################
# File Name: fetch_sess.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Tue Jul 10 10:41:06 2018
#########################################################################
#!/bin/bash

mv deci-test1b.arc deci-test1b.arc.bk
mv deci-test2b.arc deci-test2b.arc.bk

for line in `cat session-list1b`
do
	if [ ${#line}==9 ]
	then
		grep "${line} " opa-sx-180425.arc >>deci-test1b.arc
	elif [ ${#line}==10 ]
	then
		grep "${line}" opa-sx-180425.arc >>deci-test1b.arc
	else
		echo "couldn't find ${line}"
	fi
done

for line in `cat session-list2b`
do
	if [ ${#line}==9 ]
	then
		grep "${line} " opa-sx-180425.arc >>deci-test2b.arc
	elif [ ${#line}==10 ]
	then
		grep "${line}" opa-sx-180425.arc >>deci-test2b.arc
	else
		echo "couldn't find ${line}"
	fi
done
