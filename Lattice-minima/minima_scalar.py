#!/usr/bin/env python3

def random_ideals(p, N):
    import os, subprocess, threading, time, ast

    num = os.cpu_count()

    mtx = threading.Lock()
    lines = []

    def fun():
        with subprocess.Popen(
                    ['./scalar.py', str(p), str(N)],
                    stdout = subprocess.PIPE,
                    stderr = subprocess.DEVNULL,
                    env = os.environ | {'PYTHONOPTIMIZE': '1'},
                ) as proc:

            while True:
                try:
                    line = proc.stdout.readline().decode().strip()
                except Exception as e:
                    print(f'\x1b[31mthread crashed: {e}\x1b[0m')
                    break
                with mtx:
                    lines.append(line)

    threads = []
    while len(threads) < num:
        threads.append(threading.Thread(target=fun, daemon=True))
    for t in threads:
        t.start()
    
    while True:
        time.sleep(float(.1))
        with mtx:
            while lines:
                line = lines.pop(0)
                print(line, flush=True)
                yield line

def main(p, N, howmany=None):
    if howmany is None:
        howmany = 10
    results = list(zip(range(howmany), random_ideals(p, N)))

if __name__ == '__main__':
    import sys
    p = int(sys.argv[1])
    N = int(sys.argv[2])
    howmany = int(sys.argv[3]) if len(sys.argv) > 3 else None
    main(p, N, howmany)
