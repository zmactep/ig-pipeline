INSTALLATION
============
* Make sure that ../svm_data_generator is compiled! (cd ../svm_data_generator && compile.sh).
* Install protobuf compiler version 2.5.0.
* fix tools_root in src/main/resources/application.conf to point to ./train_model_example.sh dir
* mvn clean package

RUN
===
java -jar ./target/ig-backend-1.0-SNAPSHOT.jar

USAGE
=====
GET:
curl http://localhost:9999/?query=whatever
POST:
curl -F query=test -F another=testtest http://localhost:9999/
You will get job id in response. You can now pass it as a param:
curl http://localhost:9999/?query=get_result:%200
curl -F 'query=get_result: 0' http://localhost:9999/
output will be shown if it is ready
Logs are in ig-backend.log
