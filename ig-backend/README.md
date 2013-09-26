INSTALLATION
============
Make sure that ../svm_data_generator is compiled! (cd ../svm_data_generator && compile.sh). Then run:
mvn clean package

RUN
===
java -jar ./target/ig-backend-1.0-SNAPSHOT.jar

USAGE
=====
curl http://localhost:8080/?query=whatever
You will get job id in response. You can now pass it as a param:
curl http://localhost:8080/?query=get_result:%200
output will be shown if it is ready
Logs are in ig-backend.log
