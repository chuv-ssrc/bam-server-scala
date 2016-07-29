#!/usr/bin/env bash

VERSION=1.0


if [ "$#" -lt 5 ]; then
    echo ""
    echo "Usage: deploy.sh <REMOTE> <REMOTE_DIR> <PORT> <PORT_HTTPS> <SETTINGS> <SECRET>"
    echo "Params:"
    echo "  REMOTE: name of the server, e.g. varapp@varapp.vital-it.ch"
    echo "  REMOTE_DIR: path on the destination server where to copy the app archive, e.g. /home/varapp/tools/bam-server/"
    echo "  PORT: HTTP port to serve the app, e.g. 9000. Set to 'disabled' to not use HTTP."
    echo "  PORT_HTTPS: HTTPS port to serve the app, e.g. 9443. Set to 'disabled' to not use HTTPS."
    echo "  SETTINGS: settings file, e.g. /home/varapp/tools/bam-server/conf/dev.conf"
    echo "  SECRET: the secret key (overwrites what is defined in SETTINGS)."
    echo ""
    exit 1
fi;

REMOTE=$1
REMOTE_DIR=$2
PORT=$3
PORT_HTTPS=$4
SETTINGS=$5
SECRET=$6

echo ""
echo "REMOTE: "$REMOTE
echo "REMOTE_DIR: "$REMOTE_DIR
echo "PORT: "$PORT
echo "PORT_HTTPS: "$PORT_HTTPS
echo "SETTINGS: "$SETTINGS
echo "SECRET: "*
echo ""

activator dist
source="bam-server-$VERSION-SNAPSHOT"
scp target/universal/$source.zip $REMOTE:$REMOTE_DIR/
output=$(ssh -n $REMOTE "
    source ~/.bash_profile
    java -version | xargs echo 'Java version: '
    cd $REMOTE_DIR
    echo 'Kill previously running versions'
    for dir in bam-server-*/ ; do
        pid_file=\$dir/RUNNING_PID
        if [ -f \$pid_file ] ; then
            cat \$pid_file | xargs kill
            rm \$pid_file
        fi
    done;
    echo 'Remove existing sources of the same version'
    if [ -d $source ]; then
        rm -rf $source
    fi
    unzip $source.zip >/dev/null
    rm $source.zip
    cd $source
    echo 'Serving bam-server on HTTPS port $PORT...'
    nohup ./bin/bam-server -v \
        -Dplay.crypto.secret=$SECRET \
        -Dconfig.file=$SETTINGS \
        -Dhttp.port=$PORT \
        -Dhttps.port=$PORT_HTTPS \
        >/dev/null 2>&1 &

    sleep 2
    cat $REMOTE_DIR/$source/RUNNING_PID | xargs echo '... Running in process'
    echo ''
")

echo "SSH output: $output"
