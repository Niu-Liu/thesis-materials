cat=('opa2018a_ref00_GA3_0' 'opa2018a_ref00_GA3_5' 'opa2018a_ref00_GA4_0' 'opa2018a_ref00_GA4_5' 'opa2018a_ref00_GA5_0' 'opa2018a_ref00_GA5_5' 'opa2018a_ref00_GA6_0' 'opa2018a_ref00_GA6_5' 'opa2018a_ref00_GA7_0' 'opa2018a_ref00_GA7_5' 'opa2018a_ref00_GA8_0' 'opa2018a_ref00_GA8_5' 'opa2018a_ref00_GA9_0' 'opa2018a_ref00_GA9_5')

for cati in ${cat[@]};
do
	echo $cati
#	sed -n "12p" ${cati}_GaiaDR1_vsh_param.tex
	sed -n "12p" ${cati}_GaiaDR2_vsh_param.tex
#	sed -n "12p" ${cati}_icrf2_vsh_param.tex
done
