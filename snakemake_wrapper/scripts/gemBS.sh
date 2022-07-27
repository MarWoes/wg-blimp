#!/usr/bin/bash
cd $alignDir
gemBS prepare -c gemBS.conf -t gemBS_cleaned.csv -o gemBS.json
gemBS index
gemBS map
for file in *.bam; do 
	mv "$file" "$(echo "$file" | sed s/bam/unmarked.bam/)"; 
done
rm *.csi *.err *.md5
exit
