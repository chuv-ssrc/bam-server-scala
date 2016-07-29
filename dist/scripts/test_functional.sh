#!/usr/bin/env bash

printf "\n* Server is responding:\n\n"

curl localhost:9000/

printf "\n\n* Response to a dowloadRange query at localhost:9000:\n\n"

curl -I localhost:9000/downloadRange/test -H "Range: bytes=0-21793" -H "Connection: keep-alive"
echo "#Characters: " $(curl -s localhost:9000/downloadRange/test -H "Range: bytes=0-21793" -H "Connection: keep-alive" | wc -c)

printf "\n* Response to a dowloadRange query at localhost:9443:\n\n"

curl -I localhost:9443/downloadRange/test -H "Range: bytes=0-21793" -H "Connection: keep-alive"
echo "#Characters: " $(curl -s localhost:9443/downloadRange/test -H "Range: bytes=0-21793" -H "Connection: keep-alive" | wc -c)

printf '\n'