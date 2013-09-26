Instructions
============

To run this stuff just do

'''
    pip install -r requirements.txt   # To install all additional libs
    ./run_pysmtp.py                   # To run development SMTP-server (for registration activation)
    python test_server.py             # To run socket server on localhost:9999
    python manage.py syncdb           # To create sqlite3 database and superuser
    python manage.py runserver        # To run development server
'''