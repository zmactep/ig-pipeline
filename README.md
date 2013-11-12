ig-pipeline
===========

Зависимости
===========
Для тулов и базы необходимы:

1) BioPython

2) scipy/numpy

3) matplotlib

4) clustalo

5) MySQLdb для python3

6) Зависимости в ig-frontend/requirements.txt

7) GCC от 4.7 и выше. Перед сборкой нужно убедиться, что в ig-tools/ig-snooper/svm_data_generator скрипт ./rebuild_all.sh отработает без ошибок - это будет означать, что нужные библиотеки/компиляторы для C++ в системе есть.

8) protoc 2.4.1 (protobuf compiler)

9) Что-то еще, что нужно будет дописать сюда

Все проверялось только с python3

Подготовка базы
===============

1) Для MySQL выполнить скрипт ig-db/create_db.sql до строки "after syncdb"

2) Выполнить "python manage.py syncdb" и "python manage.py syncdb --database=ig" в ./ig-frontend

3) Выполнить в скрипте ig-db/create_db.sql команды после "after syncdb"


Подготовка конфигов
===================
1) В ig-config/* необходимо прописать правильные настройки подключения к базе и правильный storage_root

Установка
=========

1) В директории ./pkgbuild запустить скрипт build_all.sh

2) Там же запустить install_all.sh от root

В дальнейшем, при перезаливке сервера нужно будет лишь выполнять шаги раздела "Установка".
