#!/usr/bin/env bash

printf "\n[TEST] Server is responding to HTTP request:\n\n"

curl http://localhost:9000/

printf "\n\n[TEST] Server is responding to HTTPS request:\n\n"

curl -k https://localhost:9443/

printf "\n\n[TEST] Response to http://localhost:9000/downloadRange/test :\n\n"

curl -k -I localhost:9000/downloadRange/test -H "Range: bytes=0-21793" -H "Connection: keep-alive"
echo "#Characters: " $(curl -k -s localhost:9000/downloadRange/test -H "Range: bytes=0-21793" -H "Connection: keep-alive" | wc -c)

printf "\n[TEST] Response to https://localhost:9443/downloadRange/test :\n\n"

curl -k -I localhost:9443/downloadRange/test -H "Range: bytes=0-21793" -H "Connection: keep-alive"
echo "#Characters: " $(curl -k -s localhost:9443/downloadRange/test -H "Range: bytes=0-21793" -H "Connection: keep-alive" | wc -c)

printf '\n'