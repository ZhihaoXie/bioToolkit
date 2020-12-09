#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Version:   v1.0.0
# Author:    Zhihao Xie  \(#^o^)/
# Date:      2018.07.14
# CopyRight: Copyright (c) Zhihao Xie, All rights reserved.

"""
multi task to run from a shell cmd file,
per line is a task.
"""

import sys
import os
import re
import threading
import queue
import shutil


class MyThread(threading.Thread):
    """
    线程模型
    """
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue
        self.start()

    def run(self):
        while True:
            try:
                task = self.queue.get(block=False)
                os.system(task)  # 执行任务
                self.queue.task_done()
            except Exception as e:
                break

class MyThreadPool():
    def __init__(self, queue, size):
        self.queue = queue
        self.pool = []
        for i in range(size):
            self.pool.append(MyThread(queue))

    def joinAll(self):
        for thd in self.pool:
            if thd.isAlive():
                thd.join()


def main():
    if len(sys.argv) < 2:
        print("Note: multi task to run from a shell cmd file.")
        print("Usage: python3 %s <cmd_sh> [task_number]" % (sys.argv[0]))
        sys.exit(1)

    cmd_file = os.path.abspath(sys.argv[1])
    if len(sys.argv) >= 3:
        number = int(sys.argv[2])
        if number <= 0:
            number = 3
    else:
        number = 3
    # queue
    q = queue.Queue()
    with open(cmd_file) as fh:
        for line in fh:
            if re.search(r'^\s*$', line):
                continue
            else:
                line = line.rstrip()
                q.put(line)
    pool = MyThreadPool(queue=q, size=number)
    pool.joinAll()


if __name__ == "__main__":
    main()
