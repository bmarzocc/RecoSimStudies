#!/bin/sh
proxy=$1
cp /tmp/$proxy $PWD/
export X509_USER_PROXY=$PWD/$proxy
