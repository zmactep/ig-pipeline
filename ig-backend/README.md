INSTALLATION
============
* Make sure that ../svm_data_generator is compiled! (cd ../svm_data_generator && compile.sh).
* Install protobuf compiler version 2.4.1.
* fix tools_root in src/main/resources/application.conf to point to ./train_model_example.sh dir
* mvn clean package

RUN
===
* java -jar ./target/ig-backend-1.0-SNAPSHOT.jar
* Logs are in ig-backend.log

USAGE
=====
via HTTP:
Generate model:
* curl -H 'Accept: application/json' -X POST -d '{"task" : "generate model", "input": {"files": ["/Users/Kos/Dropbox/Biocad/ig-pipeline/data/train/VDJH_train.fasta", "/Users/Kos/Dropbox/Biocad/ig-pipeline/data/train/VDJH_train.kabat"], "params": {"mlWindowsize": "13", "algo": "random forest", "algoParams": "-l 10 -S 0", "modelName": "model1"}, "comment": "I am cool!", "group": "regions"}, "output": {"outdir": "task1/"}}' http://localhost:8080/
* curl -H 'Accept: application/json' -X POST -d '{"result_for":"0"}' http://localhost:8080/

List models:
* curl -H 'Accept: application/json' -X POST -d '{"task" : "model list", "input": {"group": "regions"}}' http://localhost:8080/
* curl -H 'Accept: application/json' -X POST -d '{"result_for":"0"}' http://localhost:8080/

Predict:
* curl -H 'Accept: application/json' -X POST -d '{"task" : "find patterns", "input": {"files": ["/Users/Kos/Dropbox/Biocad/ig-pipeline/data/train/VDJH_train.fasta", "/Users/Kos/Dropbox/Biocad/ig-pipeline/data/train/VDJH_train.kabat"], "params": {"mlWindowsize": "13", "avgWidowsize": "1", "mergeThreshold": "7", "modelPath": "task1/model1"}}, "output": {"outdir": "task1/"}}' http://localhost:8080/
* curl -H 'Accept: application/json' -X POST -d '{"result_for":"0"}' http://localhost:8080/

via browser:
* Don't forget to escape curly brackets: http://localhost:8080/?query=%7B%22result_for%22:%220%22%7D

via TCP:
* connect: telnet localhost 9999
* just paste the query: {"result_for":"0"}

CONFIGURATION FILE
==================
You can _override_ existing (but not add new) values with external config. To use external config file run:
* java -jar -Dconfig.file=/absolute/path/ig-backend.conf  ./target/ig-backend-1.0-SNAPSHOT.jar

Example of ig-backend.conf (values are the same as in application.conf, but it generally must differes):

include "application"

ig-backend {
    http_port = 8080
    tcp_port = 9999
    tools_root = /Users/Kos/Dropbox/Biocad/ig-pipeline/tools/
    storage_root = /Users/Kos/Dropbox/Biocad/ig-pipeline/storage/
    db_host = localhost
    db_port = 3306
    db_name = ig
    db_user = root
    db_password = password
}