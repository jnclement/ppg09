#!/bin/bash

if [ $# -lt 1 ]; then
    echo "need radius index"
    exit 1
fi

bash makesubrun.sh mb 0 $1
bash makesubrun.sh jet5 0 $1
bash makesubrun.sh jet12 0 $1
bash makesubrun.sh jet20 0 $1
bash makesubrun.sh jet30 0 $1
bash makesubrun.sh jet40 0 $1
bash makesubrun.sh jet50 0 $1
bash makesubrun.sh jet60 0 $1
