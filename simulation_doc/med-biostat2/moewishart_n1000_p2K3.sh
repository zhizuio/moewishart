#!/bin/bash
# The script here is to open multiple screens, and send one command to one screen
# First, excute via 'bash moewishart_n1000_p2K3.sh'
# Second, uncomment the second-half script, and copy it to command line to run, so that al opened screens will be closed.

screen -dmS test01
screen -S test01 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 1\n"; 
screen -dmS test02
screen -S test02 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 2\n"; 
screen -dmS test03
screen -S test03 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 3\n"; 
screen -dmS test04
screen -S test04 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 4\n"; 
screen -dmS test05
screen -S test05 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 5\n"; 
screen -dmS test06
screen -S test06 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 6\n"; 
screen -dmS test07
screen -S test07 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 7\n"; 
screen -dmS test08
screen -S test08 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 8\n"; 
screen -dmS test09
screen -S test09 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 9\n"; 
screen -dmS test010
screen -S test010 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 10\n"; 
screen -dmS test011
screen -S test011 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 11\n"; 
screen -dmS test012
screen -S test012 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 12\n"; 
screen -dmS test013
screen -S test013 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 13\n"; 
screen -dmS test014
screen -S test014 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 14\n"; 
screen -dmS test015
screen -S test015 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 15\n"; 
screen -dmS test016
screen -S test016 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 16\n"; 
screen -dmS test017
screen -S test017 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 17\n"; 
screen -dmS test018
screen -S test018 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 18\n"; 
screen -dmS test019
screen -S test019 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 19\n"; 
screen -dmS test020
screen -S test020 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 20\n"; 
screen -dmS test021
screen -S test021 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 21\n"; 
screen -dmS test022
screen -S test022 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 22\n"; 
screen -dmS test023
screen -S test023 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 23\n"; 
screen -dmS test024
screen -S test024 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 24\n"; 
screen -dmS test025
screen -S test025 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 25\n"; 
screen -dmS test026
screen -S test026 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 26\n"; 
screen -dmS test027
screen -S test027 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 27\n"; 
screen -dmS test028
screen -S test028 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 28\n"; 
screen -dmS test029
screen -S test029 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 29\n"; 
screen -dmS test030
screen -S test030 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 30\n"; 
screen -dmS test031
screen -S test031 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 31\n"; 
screen -dmS test032
screen -S test032 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 32\n"; 
screen -dmS test033
screen -S test033 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 33\n"; 
screen -dmS test034
screen -S test034 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 34\n"; 
screen -dmS test035
screen -S test035 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 35\n"; 
screen -dmS test036
screen -S test036 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 36\n"; 
screen -dmS test037
screen -S test037 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 37\n"; 
screen -dmS test038
screen -S test038 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 38\n"; 
screen -dmS test039
screen -S test039 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 39\n"; 
screen -dmS test040
screen -S test040 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 40\n"; 
screen -dmS test041
screen -S test041 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 41\n"; 
screen -dmS test042
screen -S test042 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 42\n"; 
screen -dmS test043
screen -S test043 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 43\n"; 
screen -dmS test044
screen -S test044 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 44\n"; 
screen -dmS test045
screen -S test045 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 45\n"; 
screen -dmS test046
screen -S test046 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 46\n"; 
screen -dmS test047
screen -S test047 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 47\n"; 
screen -dmS test048
screen -S test048 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 48\n"; 
screen -dmS test049
screen -S test049 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 49\n"; 
screen -dmS test050
screen -S test050 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 50\n"; 
screen -dmS test051
screen -S test051 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 51\n"; 
screen -dmS test052
screen -S test052 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 52\n"; 
screen -dmS test053
screen -S test053 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 53\n"; 
screen -dmS test054
screen -S test054 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 54\n"; 
screen -dmS test055
screen -S test055 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 55\n"; 
screen -dmS test056
screen -S test056 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 56\n"; 
screen -dmS test057
screen -S test057 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 57\n"; 
screen -dmS test058
screen -S test058 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 58\n"; 
screen -dmS test059
screen -S test059 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 59\n"; 
screen -dmS test060
screen -S test060 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 60\n"; 
screen -dmS test061
screen -S test061 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 61\n"; 
screen -dmS test062
screen -S test062 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 62\n"; 
screen -dmS test063
screen -S test063 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 63\n"; 
screen -dmS test064
screen -S test064 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 64\n"; 
screen -dmS test065
screen -S test065 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 65\n"; 
screen -dmS test066
screen -S test066 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 66\n"; 
screen -dmS test067
screen -S test067 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 67\n"; 
screen -dmS test068
screen -S test068 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 68\n"; 
screen -dmS test069
screen -S test069 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 69\n"; 
screen -dmS test070
screen -S test070 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 70\n"; 
screen -dmS test071
screen -S test071 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 71\n"; 
screen -dmS test072
screen -S test072 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 72\n"; 
screen -dmS test073
screen -S test073 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 73\n"; 
screen -dmS test074
screen -S test074 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 74\n"; 
screen -dmS test075
screen -S test075 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 75\n"; 
screen -dmS test076
screen -S test076 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 76\n"; 
screen -dmS test077
screen -S test077 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 77\n"; 
screen -dmS test078
screen -S test078 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 78\n"; 
screen -dmS test079
screen -S test079 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 79\n"; 
screen -dmS test080
screen -S test080 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 80\n"; 
screen -dmS test081
screen -S test081 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 81\n"; 
screen -dmS test082
screen -S test082 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 82\n"; 
screen -dmS test083
screen -S test083 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 83\n"; 
screen -dmS test084
screen -S test084 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 84\n"; 
screen -dmS test085
screen -S test085 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 85\n"; 
screen -dmS test086
screen -S test086 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 86\n"; 
screen -dmS test087
screen -S test087 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 87\n"; 
screen -dmS test088
screen -S test088 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 88\n"; 
screen -dmS test089
screen -S test089 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 89\n"; 
screen -dmS test090
screen -S test090 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 90\n"; 
screen -dmS test091
screen -S test091 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 91\n"; 
screen -dmS test092
screen -S test092 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 92\n"; 
screen -dmS test093
screen -S test093 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 93\n"; 
screen -dmS test094
screen -S test094 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 94\n"; 
screen -dmS test095
screen -S test095 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 95\n"; 
screen -dmS test096
screen -S test096 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 96\n"; 
screen -dmS test097
screen -S test097 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 97\n"; 
screen -dmS test098
screen -S test098 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 98\n"; 
screen -dmS test099
screen -S test099 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 99\n"; 
screen -dmS test0100
screen -S test0100 -p 0 -X stuff "nice -15 Rscript moewishart_n1000_p2K3.R 100\n"; 

# screen -XS test01 quit;
# screen -XS test02 quit;
# screen -XS test03 quit;
# screen -XS test04 quit;
# screen -XS test05 quit;
# screen -XS test06 quit;
# screen -XS test07 quit;
# screen -XS test08 quit;
# screen -XS test09 quit;
# screen -XS test010 quit;
# screen -XS test011 quit;
# screen -XS test012 quit;
# screen -XS test013 quit;
# screen -XS test014 quit;
# screen -XS test015 quit;
# screen -XS test016 quit;
# screen -XS test017 quit;
# screen -XS test018 quit;
# screen -XS test019 quit;
# screen -XS test020 quit;
# screen -XS test021 quit;
# screen -XS test022 quit;
# screen -XS test023 quit;
# screen -XS test024 quit;
# screen -XS test025 quit;
# screen -XS test026 quit;
# screen -XS test027 quit;
# screen -XS test028 quit;
# screen -XS test029 quit;
# screen -XS test030 quit;
# screen -XS test031 quit;
# screen -XS test032 quit;
# screen -XS test033 quit;
# screen -XS test034 quit;
# screen -XS test035 quit;
# screen -XS test036 quit;
# screen -XS test037 quit;
# screen -XS test038 quit;
# screen -XS test039 quit;
# screen -XS test040 quit;
# screen -XS test041 quit;
# screen -XS test042 quit;
# screen -XS test043 quit;
# screen -XS test044 quit;
# screen -XS test045 quit;
# screen -XS test046 quit;
# screen -XS test047 quit;
# screen -XS test048 quit;
# screen -XS test049 quit;
# screen -XS test050 quit;
# screen -XS test051 quit;
# screen -XS test052 quit;
# screen -XS test053 quit;
# screen -XS test054 quit;
# screen -XS test055 quit;
# screen -XS test056 quit;
# screen -XS test057 quit;
# screen -XS test058 quit;
# screen -XS test059 quit;
# screen -XS test060 quit;
# screen -XS test061 quit;
# screen -XS test062 quit;
# screen -XS test063 quit;
# screen -XS test064 quit;
# screen -XS test065 quit;
# screen -XS test066 quit;
# screen -XS test067 quit;
# screen -XS test068 quit;
# screen -XS test069 quit;
# screen -XS test070 quit;
# screen -XS test071 quit;
# screen -XS test072 quit;
# screen -XS test073 quit;
# screen -XS test074 quit;
# screen -XS test075 quit;
# screen -XS test076 quit;
# screen -XS test077 quit;
# screen -XS test078 quit;
# screen -XS test079 quit;
# screen -XS test080 quit;
# screen -XS test081 quit;
# screen -XS test082 quit;
# screen -XS test083 quit;
# screen -XS test084 quit;
# screen -XS test085 quit;
# screen -XS test086 quit;
# screen -XS test087 quit;
# screen -XS test088 quit;
# screen -XS test089 quit;
# screen -XS test090 quit;
# screen -XS test091 quit;
# screen -XS test092 quit;
# screen -XS test093 quit;
# screen -XS test094 quit;
# screen -XS test095 quit;
# screen -XS test096 quit;
# screen -XS test097 quit;
# screen -XS test098 quit;
# screen -XS test099 quit;
# screen -XS test0100 quit;
