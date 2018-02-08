#!/bin/sh
curl -fsLO https://raw.githubusercontent.com/scijava/scijava-scripts/master/travis-build.sh
sh travis-build.sh $encrypted_f7c1b60b55b6_key $encrypted_f7c1b60b55b6_iv
