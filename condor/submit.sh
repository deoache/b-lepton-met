#!/bin/bash

export XRD_NETWORKSTACK=IPv4
export X509_USER_PROXY=$HOME/tmp/x509up
cd MAINDIRECTORY

COMMAND

xrdcp -r condor/out/ EOSDIR