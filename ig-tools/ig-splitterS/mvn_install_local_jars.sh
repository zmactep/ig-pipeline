#!/bin/bash
mvn install:install-file -Dfile=./lib/ReadSeq.jar \
                         -DgroupId=com.github.rdpstaff \
                         -DartifactId=ReadSeq \
                         -Dversion=1.0-SNAPSHOT \
                         -Dpackaging=jar