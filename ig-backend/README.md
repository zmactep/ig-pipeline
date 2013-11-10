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
* curl -H 'Accept: application/json' -X POST -d '{
                                                  "commands":[
                                                    {
                                                      "executable": "ig-snooper/train.py",
                                                      "input": {
                                                          "params": [
                                                            {"name": "fasta", "value": "/opt/ig-pipeline/data/train/VDJH_train.fasta"},
                                                            {"name": "kabat", "value": "/opt/ig-pipeline/data/train/VDJH_train.kabat"},
                                                            {"name": "model_name", "value": "model.mdl"},
                                                            {"name": "ml_window_size", "value": "5"}
                                                          ],
                                                          "comment": "ig is really cool!",
                                                          "group": "testrun"
                                                      }
                                                    },
                                                    {
                                                       "executable": "ig-snooper/predict.py",
                                                       "input": {
                                                           "params": [
                                                             {"name": "fasta", "value": "/opt/ig-pipeline/data/train/VDJH_train.fasta"},
                                                             {"name": "model_path", "value": "/opt/ig-pipeline/storage/testrun/0/model.mdl"},
                                                             {"name": "merge_threshold", "value": "1"},
                                                             {"name": "avg_window_size", "value": "10"},
                                                             {"name": "ml_window_size", "value": "5"}
                                                           ],
                                                           "comment": "ig is cool!",
                                                           "group": "testrun"
                                                       }
                                                     }
                                                  ]
                                                 }' http://localhost:8080/
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
    tools_root = /opt/ig-pipeline/tools/
    storage_root = /opt/ig-pipeline/storage/
    db_host = localhost
    db_port = 3306
    db_name = ig
    db_user = root
    db_password = password
    workers = 3
}