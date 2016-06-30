#!/usr/bin/env bash

VERSION=1.0


if [ "$#" -lt 5 ]; then
    echo ""
    echo "Usage: deploy.sh <REMOTE> <REMOTE_DIR> <PORT> <SETTINGS> <SECRET>"
    echo "Params:"
    echo "  REMOTE: name of the server, e.g. varapp@varapp.vital-it.ch"
    echo "  REMOTE_DIR: path on the destination server where to copy the app archive, e.g. /home/varapp/tools/bam-server"
    echo "  PORT: port to serve the app, e.g. 9000"
    echo "  SETTINGS: settings file, e.g. /home/varapp/tools/bam-server/conf/dev.conf"
    echo ""
    exit 1
fi;

REMOTE=$1
REMOTE_DIR=$2
PORT=$3
SETTINGS=$4
SECRET=$5

activator dist
source=bam-server-$VERSION-SNAPSHOT
scp target/universal/$source.zip $REMOTE:$REMOTE_DIR
pid=$(ssh $REMOTE "
    cd $REMOTE_DIR
    # Kill previously running versions
    for dir in bam-server-*; do
        cat \$dir/RUNNING_PID | xargs kill
        rm \$dir/RUNNING_PID
    done;
    unzip $source
    cd $source
    ./bin/bam-server -v -Dplay.crypto.secret=$SECRET -Dconfig.file=$SETTINGS -Dhttp.port=$PORT
    cat RUNNING_PID
")

echo "PID="$pid
