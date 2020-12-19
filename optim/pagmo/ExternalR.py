from subprocess import PIPE, Popen, DETACHED_PROCESS
from json import JSONDecoder, JSONEncoder
from base64 import b64decode, b64encode


class ExternalR:

    def __init__(self, cmd, cwd, args=[]):
        self.cmd = cmd
        self.cwd = cwd
        self.args = args
        self.je = JSONEncoder(ensure_ascii=True)
        self.jd = JSONDecoder()

    def submit(self, sendData=None):
        p = Popen(args=[self.cmd] + self.args, stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=self.cwd)
        output, err = p.communicate(input=sendData)
        pResult = p.wait()
        if (pResult != 0):
            raise Exception("The process exited with an error: {}".format(err))

        # The process is expected to output a JSON.
        return self.jd.decode(b64decode(output).decode())
    
    def fitness(self, x):
        j = self.je.encode(x.tolist())
        b = b64encode(j.encode())
        return self.submit(sendData=b)


class ExternalR2:
    
    def __init__(self, cmd, cwd, args=[]):
        self.cmd = cmd
        self.cwd = cwd
        self.args = args
        self.je = JSONEncoder(ensure_ascii=True)
        self.jd = JSONDecoder()
        self.proc = None
        self.output = None
        self.err = None

    def start(self, sendData=None):
        self.proc = Popen(args=[self.cmd] + self.args, stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=self.cwd)
        self.output, self.err = self.proc.communicate(input=sendData)
        return self
    
    def start_fitness(self, x):
        j = self.je.encode(x.tolist())
        b = b64encode(j.encode())
        return self.start(sendData=b)
    
    def wait(self):
        pResult = self.proc.wait()
        if (pResult != 0):
            raise Exception("The process exited with an error: {}".format(self.err))

        # The process is expected to output a JSON.
        return self.jd.decode(b64decode(self.output).decode())
