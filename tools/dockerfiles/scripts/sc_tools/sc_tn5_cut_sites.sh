#!/bin/bash
echo "Adjusting fragments data to 1bp length Tn5 cut sites"

# 10x fragments file is saved in BED format - 0-start,
# half-open (0-based). For example, for an interval of
# length 4 with coordinates 5-9 we want to have only
# its start and end - 5-6 and 8-9 correspondingly.
# -----====-       -----=--=-     
# 0123456789       0123456789

zcat $1 | awk '{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\n"$1"\t"$3-1"\t"$3"\t"$4"\t"$5}' | sort -k1,1 -k2,2n -k3,3n | bgzip > $2
tabix -p bed $2