#########################################################################
# File Name: fetch_sess.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Tue Jul 10 10:41:06 2018
#########################################################################
#!/bin/bash

if [ -e deci-test1a.arc ]
then
	mv deci-test1a.arc deci-test1a.arc.bk
fi
if [ -e deci-test2a.arc ]
then
	mv deci-test2a.arc deci-test2a.arc.bk
fi

#for line in `cat session-list1`
for line in `cat session-list1a`
do
	if [ ${#line}==9 ]
	then
		grep "${line} " opa-sx-180425.arc >>deci-test1a.arc
	elif [ ${#line}==10 ]
	then
		grep "${line}" opa-sx-180425.arc >>deci-test1a.arc
	else
		echo "couldn't find ${line}"
	fi
done

#for line in `cat session-list2`
for line in `cat session-list2a`
do
	if [ ${#line}==9 ]
	then
		grep "${line} " opa-sx-180425.arc >>deci-test2a.arc
	elif [ ${#line}==10 ]
	then
		grep "${line}" opa-sx-180425.arc >>deci-test2a.arc
	else
		echo "couldn't find ${line}"
	fi
done
