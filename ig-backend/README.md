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
server will runs train_model_example.sh (blocking call) and outputs result.
Logs are in ig-backend.log
