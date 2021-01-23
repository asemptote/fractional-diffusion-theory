#!/bin/bash

# quick scripts

JOBRATE_PYTHON=$(cat <<END
from sys import argv
import matplotlib.pyplot as plt
alphas = []
rates = []
for line in argv[1].split('\n')[2:]:  # [2:] to start from alpha=1.15
    alphas.append(float(line.split('.')[2][:3])/100)
    rates.append(float(line.split(':')[1][:-2]))
plt.plot(alphas, rates)
plt.gca().set(xlabel=r'$\alpha$', ylabel=r'Rate (Hz)')
plt.show()
END
)

case $1 in
"jobrate")  # plot the firing rate from the job output log
    jobrate=$(paste <(ls wardakgong.o$2*) <(ls wardakgong.o$2* \
        | xargs -i sh -c "tail -n 4 {} | head -1"))
    echo "$jobrate"
    python -c "$JOBRATE_PYTHON" "$jobrate"    
    unset jobrate
    ;;
*)
    echo "Usage: $0 [jobrate|]"
    ;;
esac
