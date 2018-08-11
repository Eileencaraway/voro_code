cat dumpfile.0100000000.txt |tail -n 1000 | awk '{print $1,$2,$3,$4,$14}' > cleanfile.txt
