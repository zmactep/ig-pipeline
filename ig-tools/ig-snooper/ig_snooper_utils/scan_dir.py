import argparse
import MySQLdb
from os import listdir
from os.path import isdir, isfile, basename, join


def main():
    parser = argparse.ArgumentParser(description='scans the dir and adds new files to DB')
    parser.add_argument('--host', nargs=1, help='mysql host')
    parser.add_argument('--port', nargs=1, help='mysql port')
    parser.add_argument('--db', nargs=1, help='db name')
    parser.add_argument('--user', nargs=1, help='db user')
    parser.add_argument('--password', nargs=1, help='db password')
    parser.add_argument('--dir', nargs=1, help='directory to scan')
    args = parser.parse_args()

    scan(args.host.pop(0), int(args.port.pop(0)), args.db.pop(0), args.user.pop(0), args.password.pop(0), args.dir.pop(0))


def scan(host, port, db, user, password, dir):
    if not isdir(dir):
        print('%s is not a dir. Abort.' % dir)
        exit(1)

    try:
        db = MySQLdb.connect(host=host, port=port, user=user, passwd=password, db=db, charset='utf8')
    except Exception as e:
        print('Error connecting to db: %s' % e)
        exit(1)

    cursor = db.cursor()
    sql = """SELECT file_id, comment, path FROM igstorage_storageitem"""
    cursor.execute(sql)
    already_in_db = {path for (file_id, comment, path) in cursor.fetchall()}
    add_files(cursor, dir, already_in_db)
    db.close()


def add_files(cursor, dir, already_in_db):
    print('Scanning dir %s' % dir)
    if isfile(dir):
        print('%s is a file, not a dir. Stop scanning dir' % dir)
        return

    for path in listdir(dir):
        full_path = join(dir, path)
        if isdir(full_path):
            add_files(cursor, full_path, already_in_db)
        else:
            if path not in already_in_db:
                print('Adding file %s to db' % path)
                already_in_db.add(path)
                try:
                    cursor.execute("insert into ig.igstorage_storageitem(file_id, path) values('%s', '%s')"
                               % (basename(path), full_path))
                except Exception as e:
                    print('Error adding file %s: $s' % (path, e))
            else:
                print('File %s already in db.' % path)


if __name__ == "__main__":
    main()

