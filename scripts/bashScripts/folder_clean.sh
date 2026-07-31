echo -ne "Clean the old files in folders ... at "
date

rdataDIR="../../rdata"
if [ "$(ls -A $rdataDIR)" ]; then
    rm -r ${rdataDIR}/*
fi

outputDIR="../../output"
if [ "$(ls -A $outputDIR)" ]; then
    rm -r ${outputDIR}/*
fi

plotsDIR="../../plots"
if [ "$(ls -A $plotsDIR)" ]; then
    rm -r ${plotsDIR}/*
fi

#logsDIR="../logs"
#if [ "$(ls -A $logsDIR)" ]; then
#    rm -r ${logsDIR}/*
#fi



