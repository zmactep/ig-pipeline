__author__ = 'mactep'

import socketserver


class MyTCPHandler(socketserver.BaseRequestHandler):

    def handle(self):
        data = self.request.recv(1024).strip()
        print("Some input: %s" % data.decode())
        to_send = data
        self.request.send(to_send)

if __name__ == "__main__":
    HOST, PORT = "localhost", 9999

    server = socketserver.TCPServer((HOST, PORT), MyTCPHandler)
    server.serve_forever()
