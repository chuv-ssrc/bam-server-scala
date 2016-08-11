#!/usr/bin/env bash

if [ "$#" -lt 4 ]; then
    echo ""
    echo "Usage: restart.sh <PORT> <PORT_HTTPS> <SETTINGS> <SECRET>"
    echo "Params:"
    echo "  PORT: HTTP port to serve the app, e.g. 9000. Set to 'disabled' to not use HTTP."
    echo "  PORT_HTTPS: HTTPS port to serve the app, e.g. 9443. Set to 'disabled' to not use HTTPS."
    echo "  SETTINGS: settings file, e.g. /home/varapp/tools/bam-server/conf/dev.conf"
    echo "  SECRET: the secret key (overwrites what is defined in SETTINGS)."
    echo ""
    exit 1
fi;

PORT=$1
PORT_HTTPS=$2
SETTINGS=$3
SECRET=$4

echo ""
echo "PORT: "$PORT
echo "PORT_HTTPS: "$PORT_HTTPS
echo "SETTINGS: "$SETTINGS
echo "SECRET: "*
echo ""


echo 'Kill previously running instance'
pid_file=RUNNING_PID
if [ -f \$pid_file ] ; then
    cat \$pid_file | xargs kill
    rm \$pid_file
fi

echo 'Serving bam-server on ports '$PORT'/'$PORT_HTTPS
nohup ./bin/bam-server -v \
    -Dplay.crypto.secret=$SECRET \
    -Dconfig.file=$SETTINGS \
    -Dhttp.port=$PORT \
    -Dhttps.port=$PORT_HTTPS \
    >/dev/null 2>&1 &

sleep 2
cat $pid_file | xargs echo '... Running in process'
echo ''