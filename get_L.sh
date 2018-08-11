cat dumpfile.0100000000.txt |less|grep 'ITEM: BOX' -A 1|tail -n 1|awk '{print $2}' 
