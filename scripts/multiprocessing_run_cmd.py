#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@Filename:    multiprocessing_run_cmd.py
@Version:     V0.0.0.1
@Date:        2021/08/06 15:50:39
@Description: 多进程执行shell命令
'''

import sys
import os
import re
from multiprocessing import Pool
from loguru import logger

# 初始化logger, 作为全局变量
logger.remove()
logger.add(sys.stderr, format='{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}')


def run_os_cmd(cmd_str):
    # 单纯地执行系统命令
    if cmd_str != "":
        logger.info("begin do cmd: %s" % cmd_str)
        status = os.system(cmd_str)
        if status == 0:
            logger.info("finished do cmd: %s" % cmd_str)
            return True
        else:
            logger.info("failed do cmd: %s" % cmd_str)
            return False


def multithread_run(cmd_list, pool_size=2):
    """多进程执行命令行脚本，默认进程池大小为2
    """
    pool = Pool(pool_size)
    logger.info("begin execute command by multi process")
    for cmd in cmd_list:
        pool.apply_async(run_os_cmd, args=(cmd, ))
    pool.close()
    pool.join()
    logger.info("finished execute command by multi process")


def main():
    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
        sys.stderr.write(__doc__ + "\n")
        sys.stderr.write(f"Usage: python {sys.argv[0]} cmd.sh [threads_num]\n")
        sys.exit()
    shell_cmd = os.path.abspath(sys.argv[1])
    if len(sys.argv) >= 3:
        multi_num = int(sys.argv[2])
    else:
        multi_num = 2
    cmd_list = []
    with open(shell_cmd, 'r') as f:
        for line in f:
            if re.search(r'^\s*$|^#', line):  # 空行和#开头行跳过
                continue
            elif len(line) == 0:
                break
            cmd_list.append(line.strip())
    multithread_run(cmd_list, multi_num)


if __name__ == '__main__':
    main()
